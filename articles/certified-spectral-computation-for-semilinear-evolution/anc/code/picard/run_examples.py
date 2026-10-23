"""Generate coupled examples, independently replay them, and write tables."""

from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import platform
import subprocess
import sys
from benchmarks import coefficient_error, exact_nilpotent_solution
from record_format import make_input, normalize_input


def cases():
    return {
        "P1-constant":([[0]],[[F(1,2)]],True),
        "P2-forcing":([[0,0],[0,0]],[[F(1,8),F(-1,16)],[F(1,32),F(1,64)],
                          [F(-1,128),F(1,128)],[F(1,512),0]],True),
        "P3-nilpotent2":([[0,F(1,8)],[0,0]],[[F(1,4),F(-1,8)],[F(1,16),F(1,32)]],True),
        "P4-nilpotent4":([[F(1,8) if j == i+1 else F(0) for j in range(4)] for i in range(4)],
                          [[F((-1)**(i+ell),16*8**ell) for i in range(4)] for ell in range(3)],True),
        "P5-dense5":([[F((-1)**(i+j),40) for j in range(5)] for i in range(5)],
                      [[F((-1)**(i+ell),16*8**ell) for i in range(5)] for ell in range(5)],False),
        "P6-den96":([[F(1,16),F(-1,16*((1<<96)+1))],[F(1,32),F(-1,32)]],
                     [[F(1,3),F(-1,3*((1<<96)+1))],[F(1,24),F(1,48)]],False),
        "P7-boundary":([[F(1,16),F(-1,16)],[F(-1,32),F(3,32)]],
                        [[F(1,4),F(1,4)],[F(1,8),F(-1,8)]],False),
        "P8-highdegree":([[F(1,8)]],[[F(1,8)]]+[[F(0)] for _ in range(29)]+[[F(1,1<<32)]],False),
    }


def upward(value, digits=4):
    """An exact upward-rounded scientific display for a nonnegative rational."""
    value = F(value)
    if value == 0:
        return "$0$"
    exponent = len(str(value.numerator))-len(str(value.denominator))
    scale = F(10**exponent) if exponent >= 0 else F(1,10**(-exponent))
    while value < scale:
        scale /= 10
        exponent -= 1
    while value >= 10*scale:
        scale *= 10
        exponent += 1
    unit = scale/(10**(digits-1))
    quotient = value/unit
    rounded = (quotient.numerator+quotient.denominator-1)//quotient.denominator
    if rounded == 10**digits:
        rounded //= 10
        exponent += 1
    text = str(rounded)
    return "$"+text[0]+"."+text[1:]+r"\times10^{"+str(exponent)+"}$"


def main():
    source = Path(__file__).resolve().parent
    results = source.parents[1]/"results"/"picard"
    results.mkdir(parents=True,exist_ok=True)
    rows = []
    for name,(matrix,polynomial,exact_case) in cases().items():
        external = make_input(matrix,polynomial)
        input_path = results/(name+".input.json")
        input_path.write_text(json.dumps(external,indent=2)+"\n",encoding="utf-8")
        _, rational_matrix, rational_polynomial, _ = normalize_input(external)
        benchmark = exact_nilpotent_solution(rational_matrix,rational_polynomial) if exact_case else None
        for m in (1,6,12,20):
            stem = f"{name}-m{m}"
            path = results/(stem+".record.json")
            report_path = results/(stem+".report.json")
            command = [sys.executable,"-B",str(source/"produce.py"),"--input",str(input_path),
                       "--m",str(m),"--output",str(path)]
            production = json.loads(subprocess.run(command,check=True,text=True,capture_output=True).stdout)
            command = [sys.executable,"-B",str(source/"verify.py"),str(path),"--input",str(input_path),
                       "--expected-m",str(m),"--report",str(report_path)]
            report = json.loads(subprocess.run(command,check=True,text=True,capture_output=True).stdout)
            row = {"case":name,"m":m,"input_file":input_path.name,"record_file":path.name,
                   "report_file":report_path.name,"input_bytes":input_path.stat().st_size,
                   "record_bytes":path.stat().st_size,"producer_seconds":production["producer_seconds"],
                   "verifier_seconds":report["verifier_seconds"],
                   "record_sha256":hashlib.sha256(path.read_bytes()).hexdigest()}
            for key in ("dimension","source_degree","working_bits","degree","fractional_bits",
                        "iterations","max_strong_norm","max_exact_image_norm","max_weak_step_error",
                        "max_strong_step_difference","certified_path_error","uniform_path_error","target",
                        "terminal_centers","terminal_intervals"):
                row[key] = report[key]
            if benchmark is not None:
                record = json.loads(path.read_text(encoding="utf-8"))
                actual_bound = coefficient_error(record,benchmark)
                assert actual_bound <= F(report["certified_path_error"]["exact"])
                row["exact_polynomial_error_bound"] = str(actual_bound)
                row["exact_terminal"] = [str(sum(v[j] for v in benchmark)) for j in range(len(matrix))]
            rows.append(row)
            print(stem,"ACCEPT",f'{row["producer_seconds"]:.3f}/{row["verifier_seconds"]:.3f}s',flush=True)
    summary = {"python":platform.python_version(),"platform":platform.platform(),
               "arithmetic":"exact rational and integer; decimal bounds are descriptive only",
               "rows":rows}
    (results/"summary.json").write_text(json.dumps(summary,indent=2)+"\n",encoding="utf-8")
    lines = ["% Generated by run_examples.py; displayed bounds rounded upward.",
             r"\newcommand{\PicardResultRows}{%"]
    for row in rows:
        if row["m"] != 12:
            continue
        lines.append(" & ".join([row["case"],str(row["dimension"]),str(row["source_degree"]),
                     str(row["degree"]),upward(row["certified_path_error"]["exact"]),
                     str(row["record_bytes"]),f'{row["producer_seconds"]:.3f}/{row["verifier_seconds"]:.3f}'])+r" \\")
    lines.append("}")
    (results/"table_rows.tex").write_text("\n".join(lines)+"\n",encoding="utf-8")


if __name__ == "__main__":
    main()
