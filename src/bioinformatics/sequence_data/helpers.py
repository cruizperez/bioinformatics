def wrap_sequence_string(sequence_string: str, line_width: int) -> str:
    """
    Splits a string every given number of characters (line width) and re-joins it with a newline character.
    This effectively allows to wrap long string stretches into multiple lines.

    Args:
        sequence_string (str): String (sequence) to wrap.
        line_width (int): Width of the resulting line (number of characters per line).

    Returns:
        str: Wrapped string.
    """
    return "\n".join([sequence_string[i : i + line_width] for i in range(0, len(sequence_string), line_width)])
