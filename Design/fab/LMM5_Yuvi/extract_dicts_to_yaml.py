#!/usr/bin/env python
"""
Script to extract dictionaries from Jupyter notebooks and save them as YAML files.
"""

import os
import sys
import json
import re
import ast
from collections import defaultdict

# Try to import required packages, install if not available
try:
    import yaml
except ImportError:
    print("Installing PyYAML package...")
    import subprocess
    subprocess.check_call([sys.executable, "-m", "pip", "install", "pyyaml"])
    import yaml

class DictionaryExtractor(ast.NodeVisitor):
    """
    AST visitor to extract dictionary assignments from Python code.
    """
    def __init__(self):
        self.dictionaries = {}
        self.current_scope_vars = {}
        
    def visit_Assign(self, node):
        """Visit assignment nodes to find dictionary assignments"""
        # Only process assignments with a single target (name = value)
        if len(node.targets) == 1 and isinstance(node.targets[0], ast.Name):
            target_name = node.targets[0].id
            
            # Check if the value is a dictionary
            if isinstance(node.value, ast.Dict):
                # Extract dictionary content
                dict_content = {}
                for key, value in zip(node.value.keys, node.value.values):
                    # Handle different types of keys
                    if isinstance(key, ast.Str):
                        key_name = key.s
                    elif isinstance(key, ast.Name):
                        key_name = key.id
                    else:
                        # Skip complex keys
                        continue
                    
                    # Handle different types of values
                    if isinstance(value, ast.Num):
                        dict_content[key_name] = value.n
                    elif isinstance(value, ast.Str):
                        dict_content[key_name] = value.s
                    elif isinstance(value, ast.NameConstant):
                        dict_content[key_name] = value.value
                    elif isinstance(value, ast.Name):
                        # Try to resolve variable references
                        if value.id in self.current_scope_vars:
                            dict_content[key_name] = self.current_scope_vars[value.id]
                        else:
                            dict_content[key_name] = f"<variable: {value.id}>"
                    elif isinstance(value, ast.BinOp):
                        # Try to handle simple binary operations
                        dict_content[key_name] = f"<expression>"
                    else:
                        dict_content[key_name] = f"<complex value>"
                
                self.dictionaries[target_name] = dict_content
                self.current_scope_vars[target_name] = dict_content
            
            # Check if the value is a call to dict()
            elif isinstance(node.value, ast.Call) and isinstance(node.value.func, ast.Name) and node.value.func.id == 'dict':
                dict_content = {}
                
                # Process keyword arguments (dict(key=value))
                for keyword in node.value.keywords:
                    key_name = keyword.arg
                    value = keyword.value
                    
                    # Handle different types of values
                    if isinstance(value, ast.Num):
                        dict_content[key_name] = value.n
                    elif isinstance(value, ast.Str):
                        dict_content[key_name] = value.s
                    elif isinstance(value, ast.NameConstant):
                        dict_content[key_name] = value.value
                    elif isinstance(value, ast.Name):
                        # Try to resolve variable references
                        if value.id in self.current_scope_vars:
                            dict_content[key_name] = self.current_scope_vars[value.id]
                        else:
                            dict_content[key_name] = f"<variable: {value.id}>"
                    else:
                        dict_content[key_name] = f"<complex value>"
                
                self.dictionaries[target_name] = dict_content
                self.current_scope_vars[target_name] = dict_content
            
            # Check if the value is a copy or update of another dictionary
            elif isinstance(node.value, ast.Call) and isinstance(node.value.func, ast.Name) and node.value.func.id == 'copy':
                if len(node.value.args) == 1 and isinstance(node.value.args[0], ast.Name):
                    source_dict_name = node.value.args[0].id
                    if source_dict_name in self.dictionaries:
                        self.dictionaries[target_name] = self.dictionaries[source_dict_name].copy()
                        self.current_scope_vars[target_name] = self.dictionaries[target_name]
        
        # Continue visiting child nodes
        self.generic_visit(node)

def extract_dictionaries(notebook_path):
    """
    Extract dictionaries from a Jupyter notebook using AST parsing.
    
    Args:
        notebook_path (str): Path to the Jupyter notebook file
        
    Returns:
        dict: Dictionary of extracted dictionaries with their names as keys
    """
    print(f"Processing notebook: {notebook_path}")
    
    # Read the notebook file
    with open(notebook_path, 'r', encoding='utf-8') as f:
        notebook_content = json.load(f)
    
    # Dictionary to store extracted dictionaries
    extracted_dicts = {}
    
    # Process each cell in the notebook
    for cell in notebook_content.get('cells', []):
        if cell.get('cell_type') == 'code':
            source = ''.join(cell.get('source', []))
            
            # Try to parse the cell content using AST
            try:
                tree = ast.parse(source)
                extractor = DictionaryExtractor()
                extractor.visit(tree)
                
                # Add extracted dictionaries to the result
                for dict_name, dict_content in extractor.dictionaries.items():
                    extracted_dicts[dict_name] = dict_content
            except SyntaxError:
                # If AST parsing fails, fall back to regex-based extraction
                try:
                    # Regular expressions to match dictionary definitions
                    dict_pattern = re.compile(r'(\w+)\s*=\s*(?:dict\(|\{)(.*?)(?:\)|\})', re.DOTALL)
                    
                    # Extract dictionaries defined with dict() or {}
                    for match in dict_pattern.finditer(source):
                        dict_name = match.group(1)
                        dict_content = match.group(2)
                        
                        # Store the raw content for manual inspection
                        extracted_dicts[dict_name] = f"Raw content: {dict_content}"
                except Exception as e:
                    print(f"Error in regex fallback for cell: {e}")
    
    return extracted_dicts

def save_as_yaml(dictionaries, output_dir, notebook_name):
    """
    Save extracted dictionaries as YAML files.
    
    Args:
        dictionaries (dict): Dictionary of dictionaries to save
        output_dir (str): Directory to save YAML files
        notebook_name (str): Name of the notebook (used for file naming)
    """
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Save each dictionary as a separate YAML file
    for dict_name, dict_obj in dictionaries.items():
        # Skip dictionaries that couldn't be parsed properly
        if isinstance(dict_obj, str) and dict_obj.startswith("Could not parse"):
            print(f"Skipping {dict_name} as it couldn't be parsed properly")
            continue
        
        # Create filename
        filename = f"{notebook_name}_{dict_name}.yaml"
        filepath = os.path.join(output_dir, filename)
        
        # Save as YAML
        try:
            with open(filepath, 'w', encoding='utf-8') as f:
                yaml.dump(dict_obj, f, default_flow_style=False, sort_keys=False)
            print(f"Saved {dict_name} to {filepath}")
        except Exception as e:
            print(f"Error saving {dict_name} to YAML: {e}")

def process_notebooks(notebook_paths, output_dir="yaml_output"):
    """
    Process multiple notebooks and extract dictionaries.
    
    Args:
        notebook_paths (list): List of paths to Jupyter notebooks
        output_dir (str): Directory to save YAML files
    """
    for notebook_path in notebook_paths:
        # Extract notebook name without extension
        notebook_name = os.path.splitext(os.path.basename(notebook_path))[0]
        
        # Extract dictionaries from the notebook
        dictionaries = extract_dictionaries(notebook_path)
        
        # Save dictionaries as YAML files
        if dictionaries:
            save_as_yaml(dictionaries, output_dir, notebook_name)
            print(f"Processed {len(dictionaries)} dictionaries from {notebook_name}")
        else:
            print(f"No dictionaries found in {notebook_name}")

def main():
    # Get list of notebook files in the current directory
    notebook_files = [f for f in os.listdir('.') if f.endswith('.ipynb')]
    
    if not notebook_files:
        print("No Jupyter notebook files found in the current directory.")
        return
    
    print(f"Found {len(notebook_files)} notebook files: {', '.join(notebook_files)}")
    
    # Create output directory
    output_dir = "yaml_output"
    os.makedirs(output_dir, exist_ok=True)
    
    # Process notebooks
    process_notebooks(notebook_files, output_dir)
    
    print(f"Dictionary extraction complete. Check the '{output_dir}' directory for the YAML files.")
    
    # Create a README file with usage instructions
    readme_path = os.path.join(output_dir, "README.md")
    with open(readme_path, 'w', encoding='utf-8') as f:
        f.write("""# Extracted YAML Dictionaries

This directory contains dictionaries extracted from Jupyter notebooks and saved as YAML files.

## File Naming Convention

Files are named using the pattern: `notebook_name_dictionary_name.yaml`

## Notes

- Some complex dictionary structures may not be fully parsed
- Values marked as `<variable: name>` or `<expression>` indicate references to variables or expressions that couldn't be fully resolved
- You may need to manually edit some files to correct values

## Usage in Python

```python
import yaml

# Load a YAML file
with open('path_to_yaml_file.yaml', 'r') as f:
    data = yaml.safe_load(f)

# Now you can use the dictionary
print(data)
```
""")
    
    print(f"Created README file at {readme_path}")

if __name__ == "__main__":
    main()
