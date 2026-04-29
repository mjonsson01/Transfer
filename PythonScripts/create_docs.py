import os
import sys

# Try a robust import for clang
try:
    import clang.cindex
except ImportError:
    # If standard import fails, try to force it via the libclang path
    try:
        import libclang
        # Manually add the path to the sys.path if it's being stubborn
        site_packages = next(p for p in sys.path if 'site-packages' in p)
        sys.path.append(os.path.join(site_packages, 'clang'))
        import clang.cindex
    except Exception as e:
        print(f"Critial Error: Could not load clang modules. {e}")

def get_class_structure(file_path):
    index = clang.cindex.Index.create()
    # Adding '-fsyntax-only' makes it faster since we don't need a full compile
    tu = index.parse(file_path, args=['-x', 'c++', '-std=c++17', '-fsyntax-only'])
    
    results = {}
    
    # Debug: Check if clang actually found tokens
    found_anything = False

    for node in tu.cursor.get_children():
        # Check if the node is in the file we are actually scanning
        if node.location.file and node.location.file.name == file_path:
            if node.kind in [clang.cindex.CursorKind.CLASS_DECL, clang.cindex.CursorKind.STRUCT_DECL]:
                class_name = node.spelling or "Unnamed"
                methods = []
                
                for child in node.get_children():
                    if child.kind == clang.cindex.CursorKind.CXX_METHOD:
                        methods.append(child.spelling)
                
                results[class_name] = methods
                found_anything = True
    
    if not found_anything:
        print(f"  [!] No classes found in: {os.path.basename(file_path)}")
    return results

def build_codebase_report(root_dir, output_md):
    print(f"Scanning directory: {os.path.abspath(root_dir)}")
    with open(output_md, "w", encoding="utf-8") as f:
        f.write(f"# Codebase Structure Report\n\n")

        file_count = 0
        for root, _, files in os.walk(root_dir):
            # Skip the venv folder so we don't scan python libraries
            if 'venv' in root or '.git' in root:
                continue

            headers = [f for f in files if f.endswith(('.h', '.hpp'))]
            for h_file in headers:
                file_count += 1
                h_path = os.path.abspath(os.path.join(root, h_file))
                print(f"Checking {h_file}...")
                
                structure = get_class_structure(h_path)
                
                if structure:
                    rel_path = os.path.relpath(h_path, root_dir)
                    f.write(f"## Header: `{rel_path}`\n")
                    
                    for class_name, methods in structure.items():
                        f.write(f"### Class: `{class_name}`\n")
                        if methods:
                            for m in methods:
                                f.write(f"- `{m}()`\n")
                        else:
                            f.write("- *No methods found*\n")
                        f.write("\n")
                    f.write("---\n")
        
        if file_count == 0:
            print("No .h or .hpp files were found at all! Check your target_path.")
def build_codebase_report(root_dir, output_md):
    print(f"Scanning directory: {os.path.abspath(root_dir)}")
    with open(output_md, "w", encoding="utf-8") as f:
        f.write(f"# Codebase Structure Report\n\n")

        file_count = 0
        for root, _, files in os.walk(root_dir):
            # Skip the venv folder so we don't scan python libraries
            if 'venv' in root or '.git' in root:
                continue

            headers = [f for f in files if f.endswith(('.h', '.hpp'))]
            for h_file in headers:
                file_count += 1
                h_path = os.path.abspath(os.path.join(root, h_file))
                print(f"Checking {h_file}...")
                
                structure = get_class_structure(h_path)
                
                if structure:
                    rel_path = os.path.relpath(h_path, root_dir)
                    f.write(f"## Header: `{rel_path}`\n")
                    
                    for class_name, methods in structure.items():
                        f.write(f"### Class: `{class_name}`\n")
                        if methods:
                            for m in methods:
                                f.write(f"- `{m}()`\n")
                        else:
                            f.write("- *No methods found*\n")
                        f.write("\n")
                    f.write("---\n")
        
        if file_count == 0:
            print("No .h or .hpp files were found at all! Check your target_path.")


if __name__ == "__main__":
    target_path = "Transfer/src/" 
    build_codebase_report(target_path, "Documentation/class_map.md")
    print(f"Done! Created class_map.md")