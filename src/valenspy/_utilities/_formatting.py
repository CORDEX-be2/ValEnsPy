import re
from datetime import datetime
import pandas as pd
import calendar

def escape_regex(pattern: str) -> str:
    return re.sub(r'([.*+?^${}()|\[\]\\])', r'\\\1', pattern)

def create_named_regex(pattern_string: str) -> str:
    parts = re.split(r'(<[^>]+>)', pattern_string)  # split into literals and group names
    regex_parts = []
    for i, part in enumerate(parts):
        if part.startswith('<') and part.endswith('>'):
            group_regex = f"(?P{part}.+?)"
            regex_parts.append(group_regex)
        else:
            regex_parts.append(escape_regex(part))
    return ''.join(regex_parts)

def parse_time_period(time_string: str, format: str = None) -> tuple[datetime, datetime]:
    """
    Parses a time string and a format string into a (start_datetime, end_datetime) tuple.

    Args:
        time_string (str): The time string to parse.
        format (str): The datetime format string.

    Returns:
        Tuple[datetime, datetime]: The start and end datetime that the input represents.
    """  

    try:
        #Check if the time_string can be parsed as a pd.Period
        period = pd.Period(time_string)
        return period.start_time, period.end_time
    except ValueError:
        # If it cannot be parsed as a Period, we will try to parse using the format
        pass

    
    # Determine end time based on the most specific field present
    if format:
        # If a format is provided, but the time string is shorter than expected, we try to cut the format to match the string len
        while len(time_string) < expected_length_of_format(format):
            last_percent = format.rfind('%')
            format = format[:last_percent]
            

        start = pd.to_datetime(time_string, format=format)
        if '%S' in format:
            end = start
        elif '%M' in format:
            end = start.replace(second=59)
        elif '%H' in format:
            end = start.replace(minute=59, second=59)
        elif '%d' in format:
            end = start.replace(hour=23, minute=59, second=59)
        elif '%m' in format:
            _, last_day = calendar.monthrange(start.year, start.month)
            end = start.replace(day=last_day, hour=23, minute=59, second=59)
        elif '%Y' in format:
            end = start.replace(month=12, day=31, hour=23, minute=59, second=59)
        else:
            raise ValueError("Unsupported or insufficient time format")
        
        return start, end

    else:
        return None, None

def expected_length_of_format(format_string: str) -> int:
    """
    Calculate the expected length of a string that a given strftime format string could parse.

    Parameters:
        format_string (str): The strftime format string (e.g., "%Y-%m-%d").

    Returns:
        int: The expected length of the string.
    """
    # Map of format specifiers to their expected lengths
    format_lengths = {
        '%Y': 4,  # Year (e.g., 2023)
        '%m': 2,  # Month (e.g., 01)
        '%d': 2,  # Day (e.g., 31)
        '%H': 2,  # Hour (e.g., 23)
        '%M': 2,  # Minute (e.g., 59)
        '%S': 2,  # Second (e.g., 59)
        '%f': 6,  # Microsecond (e.g., 123456)
        '%j': 3,  # Day of the year (e.g., 001 to 366)
        '%U': 2,  # Week number of the year (Sunday as the first day of the week)
        '%W': 2,  # Week number of the year (Monday as the first day of the week)
        '%z': 5,  # UTC offset (e.g., +0100)
        '%Z': 3,  # Timezone name (e.g., UTC)
        '%a': 3,  # Abbreviated weekday name (e.g., Mon)
        '%A': 9,  # Full weekday name (e.g., Monday)
        '%b': 3,  # Abbreviated month name (e.g., Jan)
        '%B': 9,  # Full month name (e.g., January)
        '%%': 1,  # Literal '%' character
    }

    # Replace each format specifier with its expected length
    length = 0
    i = 0
    while i < len(format_string):
        if format_string[i] == '%':  # Found a format specifier
            specifier = format_string[i:i+2]
            if specifier in format_lengths:
                length += format_lengths[specifier]
                i += 2  # Skip the specifier
            else:
                raise ValueError(f"Unsupported format specifier: {specifier}")
        else:
            length += 1  # Count literal characters
            i += 1

    return length