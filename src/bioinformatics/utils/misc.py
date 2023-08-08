"""
Provides miscellaneous functions used throughout the codebase.
"""

from typing import List


def comma_separated_list_str(value: str) -> List[str]:
    """
    Converts a comma-separated string into a list.

    Args:
        value (Any): String to convert.

    Returns:
        List[Any]: _description_
    """
    if not value:
        return []
    return value.split(",")
