import os
import re
from collections import defaultdict

# Regex patterns for class/struct and functions
CLASS_PATTERN = re.compile(r'\b(class|struct)\s+(\w+)\s*[{]', re.MULTILINE)
FUNC_PATTERN = re.compile(
    r'^(?:[\w:\<\>\*\&\s]+)\s+([A-Za-z_]\w*)\s*\([^;{]*\)\s*[{;]',
    re.MULTILINE
)
CONTROL_KEYWORDS = {"if", "for", "while", "switch", "catch"}

def is_valid_function(name):
    return name not in CONTROL_KEYWORDS


def extract_signatures(file_path):
    try:
        with open(file_path, "r", errors="ignore") as f:
            content = f.read()
    except Exception as e:
        return [], [], [], str(e)

    class_defs = []
    class_methods = []
    free_functions = []

    # Find classes/structs
    for class_match in CLASS_PATTERN.finditer(content):
        class_type, class_name = class_match.groups()
        class_defs.append(f"{class_type} {class_name}")

        # Capture class body with brace matching
        start_idx = class_match.end()
        brace_count, i = 1, start_idx
        while i < len(content) and brace_count > 0:
            if content[i] == "{":
                brace_count += 1
            elif content[i] == "}":
                brace_count -= 1
            i += 1
        class_body = content[start_idx:i]

        # Extract method signatures inside class body
        for m in FUNC_PATTERN.finditer(class_body):
            signature = m.group(0).strip()
            class_methods.append(f"{class_name}::{signature}")

    # Extract top-level/free functions (not inside classes)
    for m in FUNC_PATTERN.finditer(content):
        func_name = m.group(1)
        if not is_valid_function(func_name):
            continue
        signature = m.group(0).strip()
        # if not any(signature in cm for cm in class_methods):
        #     free_functions.append(signature)
        # Use helpers to decide what to keep
        if is_function_declaration(signature, file_path) or is_function_definition(signature, file_path):
            # avoid duplication if already captured as class method
            if not any(signature in cm for cm in class_methods):
                # strip '{' for cleaner cpp output
                clean = signature.rstrip("{").strip()
                clean = signature.rstrip("{").strip()
                clean = normalize_signature(clean)
                if not any(normalize_signature(cm) == clean for cm in class_methods):
                    free_functions.append(clean) 

    free_functions = sorted(set(normalize_signature(f) for f in free_functions))
    return free_functions, class_defs, class_methods, None

def scan_directory(root_dir, extensions=(".c", ".cpp", ".h", ".hpp")):
    grouped_results = defaultdict(lambda: {
        "files": set(),
        "free_functions": set(),
        "classes": set(),
        "methods": set(),
        "errors": []
    })

    for subdir, _, files in os.walk(root_dir):
        for file in files:
            if file.endswith(extensions):
                base_name, ext = os.path.splitext(file)
                path = os.path.join(subdir, file)

                free_funcs, classes, methods, error = extract_signatures(path)

                grouped_results[base_name]["files"].add(file)
                grouped_results[base_name]["free_functions"].update(free_funcs)
                grouped_results[base_name]["classes"].update(classes)
                grouped_results[base_name]["methods"].update(methods)
                if error:
                    grouped_results[base_name]["errors"].append(f"{file}: {error}")

    return grouped_results

def is_function_declaration(line, filename):
    # in header files, prototypes end with );
    return filename.endswith(".h") and line.strip().endswith(");")

def is_function_definition(line, filename):
    # in cpp files, definitions usually end with {
    return filename.endswith(".cpp") and line.strip().endswith("{")

def normalize_signature(sig: str) -> str:
    """Remove trailing semicolons and whitespace for consistent comparison."""
    return sig.rstrip(" {;").strip()


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Extract function signatures and class/struct methods")
    parser.add_argument("directory", help="Root directory to scan")
    args = parser.parse_args()

    output = scan_directory(args.directory)
    with open("system_crawler_output.txt", "w") as out:
        for unit_name, unit in output.items():
            out.write("=" * 80 + "\n")
            out.write(f"Logical Unit: {unit_name}\n")
            out.write(f"Files: {' '.join(unit['files'])}\n\n")

            if unit["classes"]:
                out.write("Classes/Structs:\n")
                for cls in sorted(set(unit["classes"])):
                    out.write(f"  {cls}\n")
                out.write("\n")

            if unit["methods"]:
                out.write("Methods:\n")
                for method in sorted(set(unit["methods"])):
                    clean = method.replace("{", "").strip()
                    out.write(f"  {clean}\n")
                out.write("\n")

            if unit["free_functions"]:
                out.write("Free Functions:\n")
                for func in sorted(set(unit["free_functions"])):
                    clean = func.replace("{", "").strip()
                    out.write(f"  {clean}\n")
                out.write("\n")
