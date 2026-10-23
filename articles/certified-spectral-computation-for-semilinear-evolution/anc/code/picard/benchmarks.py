"""Exact polynomial solutions for zero and nilpotent matrix benchmarks."""

from fractions import Fraction


def matrix_product(left, right):
    d = len(left)
    return [[sum(left[i][k]*right[k][j] for k in range(d)) for j in range(d)]
            for i in range(d)]


def exact_nilpotent_solution(matrix, polynomial):
    """Use finite matrix powers, then check the differential equation exactly."""
    d = len(matrix)
    powers = [[[Fraction(int(i == j)) for j in range(d)] for i in range(d)]]
    for _ in range(d):
        powers.append(matrix_product(powers[-1],matrix))
    if any(v for row in powers[-1] for v in row):
        raise ValueError("matrix is not nilpotent to the asserted order")
    degree = len(polynomial)+d-2
    solution = [[Fraction(0)]*d for _ in range(degree+1)]
    for j, source in enumerate(polynomial):
        denominator = 1
        for k in range(d):
            if k:
                denominator *= j+k
            for i in range(d):
                solution[j+k][i] += sum(powers[k][i][a]*source[a] for a in range(d))/denominator
    assert solution[0] == polynomial[0]
    for ell in range(degree+1):
        derivative = [(ell+1)*solution[ell+1][i] if ell < degree else Fraction(0)
                      for i in range(d)]
        force = [(ell+1)*polynomial[ell+1][i] if ell+1 < len(polynomial) else Fraction(0)
                 for i in range(d)]
        for i in range(d):
            assert derivative[i] == sum(matrix[i][a]*solution[ell][a] for a in range(d))+force[i]
    return solution


def coefficient_error(record, exact_coefficients):
    table = record["iterates"][-1]
    d = record["dimension"]
    denominator = 1 << record["fractional_bits"]
    error = Fraction(0)
    for ell in range(max(len(table),len(exact_coefficients))):
        differences = [abs((Fraction(int(table[ell][i]),denominator) if ell < len(table) else 0)
                           -(exact_coefficients[ell][i] if ell < len(exact_coefficients) else 0))
                       for i in range(d)]
        error += max(differences)
    return error
