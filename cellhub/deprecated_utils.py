import math


DEFAULT_MEM = 4  # Default memory: 4G.


def parse2int(s, *, positive_only=False):
    """
    Strip and parse a string to int.

    Accepts optional '+' or '-' sign and decimal digits.
    Rejects floats, scientific notation, and non-digit characters.

    Parameters
    ----------
    s : Any
        Value to parse.
    positive_only : bool, default False
        If True, reject negative or zero values.

    Returns
    -------
    int
        Parsed integer value.

    Raises
    ------
    ValueError
        If parsing fails or violates positive-only constraint.
    """
    v = str(s).strip()
    if not v:
        raise ValueError("Empty string cannot be parsed to int.")

    sign = v[0] if v[0] in "+-" else ""
    digits = v[1:] if sign else v

    if not digits.isdigit():
        raise ValueError(f"Cannot parse `{s}` to int.")

    num = int(v)

    if positive_only and num <= 0:
        raise ValueError(f"Value `{num}` is not positive.")

    return num


def parse2float(s, *, finite_only=False):
    """
    Strip and parse a string to float.

    Parameters
    ----------
    s : Any
        Value to parse.
    finite_only : bool, default False
        If True, reject NaN or Inf values.

    Returns
    -------
    float
        Parsed float value.

    Raises
    ------
    ValueError
        If parsing fails or violates finite-only constraint
        (when finite_only=True and is NaN/Inf).
    """
    v = str(s).strip()
    if not v:
        raise ValueError("Empty string cannot be parsed to float.")

    try:
        f = float(v)
    except (ValueError, TypeError, AttributeError):
        raise ValueError(f"Cannot parse `{s}` to float.")

    if finite_only and not math.isfinite(f):
        raise ValueError(f"Value `{s}` is not finite (got {f}).")

    return f


def str2list(s):
    """
    Convert comma-separated string to list of strings.
    - Trims whitespace
    - Filters out empty entries
    - Accepts list input and returns as-is
    - Converts None to empty list
    """
    if s is None:
        return []
    if isinstance(s, list):
        return s
    return [x.strip() for x in s.strip().split(",") if x.strip()]


def parse_mem(memory):
    """
    Return an integer that represents the amount of memory
    needed by the task in gigabytes.
    """

    if memory in [None, False]:
        return DEFAULT_MEM

    if isinstance(memory, (int, float)):
        return int(round(memory))

    if isinstance(memory, str):
        memory = memory.strip().lower()
        if memory in ["", "none", "false", "null", "default"]:
            return DEFAULT_MEM
        try:
            if memory.endswith("g"):
                return int(round(float(memory[:-1].strip())))
            elif memory.endswith("m"):
                return int(round(float(memory[:-1].strip()) / 1000))
            elif memory.replace(".", "", 1).isdigit():
                return int(round(float(memory)))
        except ValueError:
            pass

    raise ValueError(
        f"Memory request `{memory}` not recognised.\n"
        "Please specify the memory required in gigabytes (G) or megabytes (M), "
        'e.g. "4G" or "4000M". If a unit is not specified, G will be assumed.'
    )


# if __name__ == "__main__":
#     assert parse2int("42") == 42, f"parsed result: {parse2int('42')}"
#     assert parse2int("  007 ") == 7, f"parsed result: {parse2int('  007 ')}"
#     assert parse2int(" -5  ") == -5, f"parsed result: {parse2int(' -5  ')}"
#     try:
#         assert parse2int(" -5  ", positive_only=True)
#     except ValueError:
#         pass

#     assert parse2float("   3.14 ") == 3.14, f"parsed result: {parse2float('   3.14 ')}"
#     assert parse2float("1e-3 ") == 0.001, f"parsed result: {parse2float('1e-3 ')}"
#     try:
#         parse2float("   inf ", finite_only=True)
#     except ValueError:
#         pass

#     assert parse_mem(5) == 5, f"parsed result: {parse_mem(5)}"
#     assert parse_mem(5.2) == 5, f"parsed result: {parse_mem(5.2)}"
#     assert parse_mem(None) == 4, f"parsed result: {parse_mem(None)}"
#     assert parse_mem(False) == 4, f"parsed result: {parse_mem(False)}"
#     assert parse_mem(" falsE") == 4, f"parsed result: {parse_mem(' falsE')}"
#     assert parse_mem("nOne ") == 4, f"parsed result: {parse_mem('nOne ')}"
#     assert parse_mem(" defaulT ") == 4, f"parsed result: {parse_mem(' defaulT ')}"
#     assert parse_mem("5.3G") == 5, f"parsed result: {parse_mem('5.3G')}"
#     assert parse_mem("5300 m") == 5, f"parsed result: {parse_mem(' 5300 m ')}"
#     assert parse_mem(" 4.6 g ") == 5, f"parsed result: {parse_mem(' 4.6 g ')}"
