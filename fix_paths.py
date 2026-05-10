#!/usr/bin/env python3
"""
Fix relative paths in .qmd files for resources in docs/ subfolder
Converts:
  Figure/ -> ../Figure/
  Data/ -> ../Data/
"""
import re
from pathlib import Path

docs_dir = Path('docs')

# Patterns to replace
replacements = [
    (r'src="Figure/', r'src="../Figure/'),
    (r"src='Figure/", r"src='../Figure/"),
    (r'src="Data/', r'src="../Data/'),
    (r"src='Data/", r"src='../Data/"),
    (r'href="Figure/', r'href="../Figure/'),
    (r"href='Figure/", r"href='../Figure/"),
    (r'href="Data/', r'href="../Data/'),
    (r"href='Data/", r"href='../Data/"),
]

# Process all .qmd files in docs/
for qmd_file in docs_dir.glob('*.qmd'):
    print(f"Processing {qmd_file.name}...", end=" ")

    with open(qmd_file, 'r', encoding='utf-8') as f:
        content = f.read()

    original_content = content

    # Apply all replacements
    for pattern, replacement in replacements:
        content = re.sub(pattern, replacement, content)

    # Write back if changed
    if content != original_content:
        with open(qmd_file, 'w', encoding='utf-8') as f:
            f.write(content)
        print("✓ Fixed")
    else:
        print("✓ No changes needed")

print("\n✅ Path fixing complete!")
