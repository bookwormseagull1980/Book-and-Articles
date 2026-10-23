"""Endpoint-layer definitions used by the methods section."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class ShellLayer:
    level: int
    name: str
    epoch: str
    endpoint_law: str
    status: str
    manuscript_role: str


SHELL_LAYERS = (
    ShellLayer(1, "particle", "early high-energy sector", "spin-2 spectral pole criterion", "recorded", "context"),
    ShellLayer(2, "nuclear", "BBN", "reaction rate equals expansion rate", "recorded", "context"),
    ShellLayer(3, "atomic", "recombination", "recombination rate equals expansion rate", "recorded", "context"),
    ShellLayer(4, "stellar", "first-star formation", "Jeans entry and cooling entry", "defined", "structure layer"),
    ShellLayer(5, "galactic", "galaxy formation", "Delta^2(k_NL,z)=C_NL", "defined", "structure layer"),
    ShellLayer(6, "cluster-web", "late large-scale structure", "Delta^2(k_NL,z)=C_NL", "defined", "structure layer"),
    ShellLayer(7, "lambda", "future IR endpoint", "IR fixed point", "recorded", "endpoint"),
)


def shell_layers() -> tuple[ShellLayer, ...]:
    return SHELL_LAYERS
