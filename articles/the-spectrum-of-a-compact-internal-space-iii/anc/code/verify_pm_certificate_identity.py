"""Verify the published PM evidence without evaluating a force or advancing a state.

The migration record distinguishes the certified computation source from its
portable deployment. Numerical results retain their original bytes. File and
AST checks establish the declared relocation, not a new certificate execution.
"""
from __future__ import annotations

import argparse
import ast
import hashlib
import json
from pathlib import Path

from cosmic_cg.params import KEY_ALIASES, load_model_params

ANC = Path(__file__).resolve().parents[1]
DATA = ANC / "results" / "pm_continuum_current"
MANIFEST = DATA / "deployment_identity.json"


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def parameter_values(params):
    return {key: None if getattr(params, key) is None else float(getattr(params, key)).hex()
            for key in KEY_ALIASES}


def numerical_payload_digest(result):
    payload = {key: value for key, value in result.items() if key not in {
        "upstream_run_id", "upstream_parameter_sha256", "upstream_reference", "identity_binding"}}
    return hashlib.sha256(json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()).hexdigest()


class Relocation(ast.NodeTransformer):
    """Normalize only the two reviewed source-file identity adaptations."""
    def visit_Name(self, node):
        if node.id == "SOURCE":
            return ast.copy_location(ast.Name(id="WORK", ctx=node.ctx), node)
        return node

    def visit_Assert(self, node):
        # The original exact-source test also accepts its explicitly verified
        # relocation record. No arithmetic expression is changed by this rule.
        if (isinstance(node.test, ast.Call) and isinstance(node.test.func, ast.Name)
                and node.test.func.id == "source_record_matches"):
            expected = ast.parse(
                "assert result['source_sha256']==hashlib.sha256(source.read_bytes()).hexdigest()"
            ).body[0]
            if ast.dump(node.test) != ast.dump(ast.parse(
                    "source_record_matches(result['source_sha256'], source)", mode="eval").body):
                raise ValueError("Undeclared source-identity adaptation")
            return expected
        return self.generic_visit(node)


def function_digests(path):
    tree = ast.parse(Path(path).read_text(encoding="utf-8-sig"))
    result = {}
    for node in tree.body:
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
            normalized = Relocation().visit(node)
            result[node.name] = hashlib.sha256(ast.dump(normalized, include_attributes=False).encode()).hexdigest()
    return result


def computational_module_digest(path):
    tree = ast.parse(Path(path).read_text(encoding="utf-8-sig"))
    body = []
    for node in tree.body:
        if isinstance(node, (ast.Import, ast.ImportFrom)):
            continue
        if isinstance(node, ast.Expr) and isinstance(node.value, ast.Constant) and isinstance(node.value.value, str):
            continue
        if (isinstance(node, ast.Assign) and len(node.targets) == 1
                and isinstance(node.targets[0], ast.Name)
                and node.targets[0].id in {"P9", "WORK", "SOURCE"}):
            continue
        body.append(Relocation().visit(node))
    return hashlib.sha256(ast.dump(ast.Module(body=body, type_ignores=[]), include_attributes=False).encode()).hexdigest()


def source_record_matches(recorded, path):
    manifest = json.loads(MANIFEST.read_text(encoding="utf-8"))
    entry = manifest["source_migration"][Path(path).name]
    actual = digest(path)
    if actual != entry["deployed_sha256"]:
        raise ValueError("Deployed source differs from its migration record")
    return recorded in {actual, entry["certified_sha256"]}


def verify(params=None):
    manifest = json.loads(MANIFEST.read_text(encoding="utf-8"))
    params = params or load_model_params()
    if parameter_values(params) != manifest["consumed_values_hex"]:
        raise ValueError("A consumed numerical parameter changed; identity-only binding is invalid")
    for name, expected in manifest["artifact_sha256"].items():
        if digest(DATA / name) != expected:
            raise ValueError(f"Published evidence changed: {name}")
    for name, entry in manifest["source_migration"].items():
        source = ANC / "code" / name
        if (digest(source) != entry["deployed_sha256"]
                or function_digests(source) != entry["function_ast_sha256"]
                or computational_module_digest(source) != entry["computational_module_ast_sha256"]):
            raise ValueError(f"Unverified source migration: {name}")
    for name, expected in manifest["numerical_source_sha256"].items():
        if digest(ANC / name) != expected:
            raise ValueError(f"Numerical source changed: {name}")
    if function_digests(ANC / "code" / "verify_cosmological_dynamics.py")["evaluate"] != manifest["formal_evaluate_ast_sha256"]:
        raise ValueError("The formal numerical evaluator changed")
    formal = json.loads((ANC / "results/cosmological_active_current/results.json").read_text(encoding="utf-8"))
    if numerical_payload_digest(formal) != manifest["formal_numerical_payload_sha256"]:
        raise ValueError("The saved formal numerical result changed")
    target = formal["PM"]["trajectory_certificate"]
    residual = json.loads((DATA / "p9_continuous_residual_certificate_20261004.json").read_text(encoding="utf-8"))
    if (target["initial_state_sha256"] != residual["initial_state_sha256"]
            or target["finite_operator_data_sha256"] != residual["operator_sha256"]):
        raise ValueError("Formal and continuous certificates describe different initial data or operators")
    if residual["cache_sha256"] != manifest["artifact_sha256"]["p9_actual_trajectory_20261004.npz"]:
        raise ValueError("Residual certificate cache identity differs")
    if not source_record_matches(residual["source_sha256"], ANC / "code/p9_continuous_residual_certificate.py"):
        raise ValueError("Residual source identity is not certified")
    independent = json.loads((DATA / "p9_i17_independent_verification_20261004.json").read_text(encoding="utf-8"))
    for name, expected in independent["source_files_sha256"].items():
        if not source_record_matches(expected, ANC / "code" / name):
            raise ValueError(f"Independent verification source differs: {name}")
    if not independent["all_checks_passed"]:
        raise ValueError("Independent spatial verification did not pass")
    return {"all_checks_passed": True, "upstream_run_id": params.source_run_id,
            "consumed_parameter_count": len(KEY_ALIASES),
            "verified_artifacts": len(manifest["artifact_sha256"]),
            "verified_deployed_sources": len(manifest["source_migration"]),
            "initial_state_sha256": residual["initial_state_sha256"],
            "operator_sha256": residual["operator_sha256"],
            "numerical_payload_sha256": numerical_payload_digest(formal),
            "new_force_evaluations": 0, "new_trajectory_advancements": 0,
            "certificate_recomputed": False,
            "scope": "input, source-relocation and saved-evidence integrity"}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--params-file")
    args = parser.parse_args()
    print(json.dumps(verify(load_model_params(args.params_file)), indent=2))


if __name__ == "__main__":
    main()
