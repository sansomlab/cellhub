DEFAULT_MEM = 4  # Default memory: 4G.


def parse2int(s):
    """
    Strip and parse a string to int.
    Only accepts strings that are pure digits, like "123", " 007 ".
    Rejects floats, scientific notation, and non-digit characters.
    """
    v = s.strip()
    if v.isdigit():
        return int(v)
    raise ValueError(f"Cannot parse `{s}` to int.")


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


if __name__ == "__main__":
    assert parse2int("42") == 42, f"parsed result: {parse2int('42')}"
    assert parse2int("  007 ") == 7, f"parsed result: {parse2int('  007 ')}"

    assert parse_mem(5) == 5, f"parsed result: {parse_mem(5)}"
    assert parse_mem(5.2) == 5, f"parsed result: {parse_mem(5.2)}"
    assert parse_mem(None) == 4, f"parsed result: {parse_mem(None)}"
    assert parse_mem(False) == 4, f"parsed result: {parse_mem(False)}"
    assert parse_mem(" falsE") == 4, f"parsed result: {parse_mem(' falsE')}"
    assert parse_mem("nOne ") == 4, f"parsed result: {parse_mem('nOne ')}"
    assert parse_mem(" defaulT ") == 4, f"parsed result: {parse_mem(' defaulT ')}"
    assert parse_mem("5.3G") == 5, f"parsed result: {parse_mem('5.3G')}"
    assert parse_mem("5300 m") == 5, f"parsed result: {parse_mem(' 5300 m ')}"
    assert parse_mem(" 4.6 g ") == 5, f"parsed result: {parse_mem(' 4.6 g ')}"
