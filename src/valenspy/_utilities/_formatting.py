import re

def escape_regex(pattern: str) -> str:
    return re.sub(r'([.*+?^${}()|\[\]\\])', r'\\\1', pattern)

def create_named_regex(pattern_string: str) -> str:
    parts = re.split(r'(<[^>]+>)', pattern_string)  # split into literals and group names
    regex_parts = []
    for i, part in enumerate(parts):
        if part.startswith('<') and part.endswith('>'):
            group_name = part[1:-1]
            # Look ahead to infer delimiters
            next_literal = parts[i + 1] if i + 1 < len(parts) else ''
            # Collect characters to exclude based on surrounding literals
            exclude_chars = set()
            if next_literal:
                m = re.match(r'([^\w<])', next_literal)
                if m:
                    exclude_chars.add(m.group(1))
            if not exclude_chars:
                group_regex = f"(?P<{group_name}>.+?)"
            else:
                exclude = ''.join(re.escape(c) for c in exclude_chars)
                group_regex = f"(?P<{group_name}>[^{exclude}]+)"
            regex_parts.append(group_regex)
        else:
            regex_parts.append(escape_regex(part))
    return ''.join(regex_parts)