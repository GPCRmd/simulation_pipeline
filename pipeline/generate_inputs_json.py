import json

def ask_question(prompt, default=None, is_bool=False):
    """
    Helper function to ask a question and get user input.
    If the user presses Enter without input, the default value is used.
    """
    if default is not None:
        prompt += f" (default: {default})"
    prompt += ": "

    user_input = input(prompt).strip()

    if is_bool:
        if user_input.lower() in ["yes", "y", "true", "t", "1"]:
            return True
        elif user_input.lower() in ["no", "n", "false", "f", "0"]:
            return False
        elif default is not None:
            return default
        else:
            return False

    return user_input if user_input else default


def create_entry(entries):
    """
    Create a single entry for the JSON file by asking the user questions.
    """
    print("\n--- Creating a new entry ---")
    name = ask_question("Enter the name of the system", default="system_name")
    pdbfile = ask_question("Enter the path to the PDB file", default="path/to/pdbfile.pdb")
    modres = ask_question("Enter modified residues (comma-separated, or leave blank for none)", default="")
    modres = [res.strip() for res in modres.split(",") if res.strip()]  # Convert to a list

    # Ligands
    ligands = []
    add_ligand = ask_question("Do you want to add ligands? (yes/no)", default="no", is_bool=True)
    while add_ligand:
        print("\n--- Adding a ligand ---")
        resname = ask_question("Enter the residue name of the ligand", default="LIG")
        name = ask_question("Enter the name of the ligand", default="Ligand")
        covalently_bound = ask_question("Is the ligand covalently bound? (yes/no)", default="no", is_bool=True)
        inchikey = ask_question("Enter the InChIKey of the ligand", default="")
        ligands.append({
            "resname": resname,
            "name": name,
            "covalently_bound": covalently_bound,
            "inchikey": inchikey
        })
        add_ligand = ask_question("Do you want to add another ligand? (yes/no)", default="no", is_bool=True)

    apo = ask_question("Do you want to simulate an apo version? (yes/no)", default="no", is_bool=True)
    prot_chain = ask_question("Enter the chain ID of the main protein", default="A")
    pdbcode = ask_question("Enter the closest-ressembling PDB code (or leave blank for none)", default="")
    curated = ask_question("Has the system been properly curated? (yes/no)", default="no", is_bool=True)
    sod2x50 = ask_question("Does the system require sodium near 2x50? (yes/no)", default="no", is_bool=True)
    isgpcr = ask_question("Is the system a GPCR? (yes/no)", default="yes", is_bool=True)

    # Return the entry as a dictionary
    entry = {
        "name": name,
        "pdbfile": pdbfile,
        "modres": modres,
        "ligands": ligands,
        "apo": apo,
        "prot_chain": prot_chain,
        "pdbcode": pdbcode,
        "curated": curated,
        "sod2x50": sod2x50,
        "isgpcr": isgpcr
    }
    
    entries.append(entry)
    
    # If apo is True, create an additional entry with apo=False
    if apo:
        apo_entry = entry.copy()
        apo_entry["apo"] = False
        entries.append(apo_entry)

def main():
    print("Welcome to the JSON generator for the simulation pipeline input!")

    # Ask if the user wants a simple entry or multiple entries
    multiple_entries = ask_question("Do you want to create multiple entries? (yes/no)", default="no", is_bool=True)

    entries = []

    if multiple_entries:
        # Create multiple entries
        add_entry = True
        while add_entry:
            create_entry(entries)
            add_entry = ask_question("\nDo you want to add another entry? (yes/no)", default="no", is_bool=True)
    else:
        # Create a single entry
        create_entry(entries)


    # Save to JSON file
    output_dir = ask_question("\nEnter the directory to save the output JSON file", default="../simulation_output/input_jsons")
    output_file = ask_question("Enter the name of the output JSON file", default="inputs.json")
    full_path = f"{output_dir.rstrip('/')}/{output_file}"
    with open(full_path, "w") as f:
        json.dump(entries, f, indent=4)
    print(f"\nJSON file saved as {full_path}")

if __name__ == "__main__":
    main()