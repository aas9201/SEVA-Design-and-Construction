from plasmid import Part, Module, Plasmid
import tkinter as tk
from tkinter import ttk, messagebox
import re
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio import SeqIO
import os 
import csv 
import numpy as np


def sort_sites(plasmid_sequence):
        """Sorts the imported plasmid sequence and does some structual validity checks"""
        predefined_sites = {
            "OriT": "CTTTTCCGCTGCATAACCCTGCTTCGGGGTCATTATAGCGATTTTTTCGGTATATCCATCCTTTTTCGCACGATATACAGGATTTTGCCAAAGGGTTCGTGTAGACTTTCCTTGGTGTATCCAACGGCGTCAGCCGGGCAGGATAGGTGAAGTAGGCCCACCCGCGAGCGGGTGTTCCTTCTTCACTGTCCCTTATTCGCACCTGGCGGTGCTCAACGGGAATCCTGCTCTGCGAGGCTGGCCGTA",
            "T0": "CTTGGACTCCTGTTGATAGATCCAGTAATGACCTCAGAACTCCATCTGGATTTGTTCAGAACGCTCGGTTGCCGCCGGGCGTTTTTTATTGGTGAGAATCCAG",
            "T1": "CAGCTGTCTAGGGCGGCGGATTTGTCCTACTCAGGAGAGCGTTCACCGACAAACAACAGATAAAACGAAAGGCCCAGTCTTTCGACTGAGCCTTTCGTTTTATTTGATGCCT",
            "PshAI": "GACNNNNGTC",
            "SwaI": "ATTTAAAT",
            "AscI": "GGCGCGCC",
            "FseI": "GGCCGGCC",
            "PacI": "TTAATTAA",
            "SpeI": "ACTAGT",
            "SanDI": "GGGTCCC",
        }

        reference_point = plasmid_sequence.find(
            "TTAATTAA"
        )  # This just makes PacI as the starting reference point of every plasmid sequence due to it being circular
        if reference_point == -1:
            raise ValueError("PacI RE site not found, please enter valid SEVA sequence")
        ordered_sequence = (
            plasmid_sequence[reference_point:] + plasmid_sequence[:reference_point]
        )

        correct_order = [
            "PacI",
            "SpeI",
            "T0",
            "SanDI",
            "SwaI",
            "PshAI",
            "OriT",
            "FseI",
            "AscI",
            "T1",
        ]
        site_positions = {}
        first_occurance = []

        for (
            site_name,
            site_sequence,
        ) in (
            predefined_sites.items()
        ):  # Subs in the "NNN" in the PshAI as this RE site could be any combination
            pattern = re.sub("N", "[ACGT]", site_sequence)
            positions = []

            for i in range(len(ordered_sequence)):
                if re.match(
                    pattern, ordered_sequence[i:]
                ):  # Checks if these sites are within the plasmid sequence and what positions they are in
                    positions.append(i + 1)

            if positions:
                site_positions[site_name] = positions
                first_occurance.append(
                    (site_name, positions[0])
                )  # Creates a new dictionary to store site and position

            first_occurance.sort(
                key=lambda x: x[1]
            )  # Sorts them by where they first appeared, as an invalid plasmid may contain more of one of these sites which is invalid

        sorted_sites = {site: site_positions[site] for site, _ in first_occurance}

        for i in sorted_sites:
            if len(sorted_sites[i]) > 1:
                raise ValueError(
                    f"The plasmid contains more than one: {i}"
                )  # Checks if there is more than one each non-variable site which is invalid
        
        plasmid_order = list(sorted_sites.keys())
        if (
            plasmid_order != correct_order
        ):  # Checks to see if the plasmid has correct order/structure
            raise ValueError("Please make sure plasmid is of valid SEVA Format")

        return sorted_sites
    


def get_module(plasmid_sequence, name, ID, module_type, description=""):
    """Creates a new module object by slicing the plasmid sequence depending on which module_type you want"""
    sorted_sites = sort_sites(
        plasmid_sequence
    )  # Sorts the plasmid so that it can be sliced easily

    if (
        module_type == "cargo"
    ):  # Slices the plasmid between restriction enzyme site to get Cargo module
        cargo_sequence = plasmid_sequence[
            (sorted_sites["PacI"][0] + 7) : sorted_sites["SpeI"][0] - 1
        ]
        cargo_part = Part(f"{name}_part", f"{ID}_part", cargo_sequence, "cargo")
        cargo_module = Module(name, ID, "cargo", description, [cargo_part])

        return cargo_module

    if (
        module_type == "marker"
    ):  # Slices the plasmid between restriction enzyme site to get Marker module
        marker_sequence = plasmid_sequence[
            (sorted_sites["SwaI"][0] + 7) : sorted_sites["PshAI"][0] - 1
        ]
        marker_part = Part(f"{name}_part", f"{ID}_part", marker_sequence, "abr")

        SwaI = Part("SwaI", "swai", "ATTTAAAT", "re")
        pshai = plasmid_sequence[
            (sorted_sites["PshAI"][0] - 1) : sorted_sites["PshAI"][0] + 9
        ]  # As the PshAI contains NNN this needs to be created dpending on imported sequence
        PshAI = Part("PshAI", "pshai", pshai, "re")
        marker_module = Module(name, ID, "marker", description, structure= [SwaI,marker_part,PshAI])
        return marker_module

    if (
        module_type == "replication"
    ):  # Slices the plasmid between restriction enzyme site to get Replication module
        replication_sequence = plasmid_sequence[
            (sorted_sites["FseI"][0] + 7) : sorted_sites["AscI"][0] - 1
        ]
        replication_part = Part(f"{name}_part", f"{ID}_part", replication_sequence, "ori")
        replication_module = Module(name, ID, "replication", description, [replication_part])

        return replication_module


def get_scars(
    plasmid_sequence,
):  # Slices plasmid between certain sites to check for any potential scars as each SEVA plasmid has certain scars, some don't
    sorted_sites = sort_sites(plasmid_sequence)
    scars = {}
    check_scars = {
        "Scar after T1": (plasmid_sequence[0 : sorted_sites["PacI"][0] - 1]),
        "Scar before T0": (
            plasmid_sequence[(sorted_sites["SpeI"][0] + 5) : sorted_sites["T0"][0] - 1]
        ),
        "Scar after T0": (
            plasmid_sequence[
                (sorted_sites["SanDI"][0] + 6) : sorted_sites["SwaI"][0] - 1
            ]
        ),
        "Scar before OriT": (
            plasmid_sequence[
                (sorted_sites["PshAI"][0] + 9) : sorted_sites["OriT"][0] - 1
            ]
        ),
        "Scar after OriT": (
            plasmid_sequence[
                (sorted_sites["OriT"][0] + 245) : sorted_sites["FseI"][0] - 1
            ]
        ),
        "Scar before T1": (
            plasmid_sequence[
                (sorted_sites["AscI"][0] + 7) : sorted_sites["SwaI"][0] - 1
            ]
        ),
    }

    for key, value in check_scars.items():
        scars[key] = value
    print(scars)
    return scars



def import_plasmid(sequence):

    # Create the main window
    check_scar = get_scars(sequence)
    data = None
    scar_data = {}
    root = tk.Tk()
    root.title("Plasmid Module Details")

    # Main frame for layout
    main_frame = ttk.Frame(root, padding="10")
    main_frame.pack(fill='both', expand=True)

    plasmid_frame = ttk.Frame(main_frame, padding="10")
    plasmid_frame.grid(row=0, column=0, sticky='nswe')

    # Frame for optional scar details
    scar_frame = ttk.Frame(main_frame, padding="10")
    scar_frame.grid(row=0, column=1, sticky='nswe', padx=(10, 0))



    # Helper function to create labeled entries or text areas
    def create_input(container, label, is_text=False):
        ttk.Label(container, text=label).grid(columnspan=2, sticky='w')
        if is_text:
            input_widget = tk.Text(container, height=3, width=40)
        else:
            input_widget = ttk.Entry(container, width=40)
        input_widget.grid(columnspan=2, sticky='ew')
        return input_widget

    # Creating input fields
    plasmid_name = create_input(plasmid_frame, "Plasmid Name (Required)")
    plasmid_id = create_input(plasmid_frame, "Plasmid ID (Required)")
    plasmid_description = create_input(plasmid_frame, "Plasmid Description (Not Required)", is_text=True)

    cargo_name = create_input(plasmid_frame, "Cargo Name (Required)")
    cargo_id = create_input(plasmid_frame, "Cargo ID (Required)")
    cargo_description = create_input(plasmid_frame, "Cargo Description (Not Required)", is_text=True)

    marker_name = create_input(plasmid_frame, "Marker Name (Required)")
    marker_id = create_input(plasmid_frame, "Marker ID (Required)")
    marker_description = create_input(plasmid_frame, "Marker Description (Not Required)", is_text=True)

    origin_name = create_input(plasmid_frame, "Origin Name (Required)")
    origin_id = create_input(plasmid_frame, "Origin ID (Required)")
    origin_description = create_input(plasmid_frame, "Origin Description (Not Required)", is_text=True)
    required_fields = ["plasmid_name", "plasmid_id", "cargo_name", "cargo_id", "marker_name", "marker_id", "origin_name", "origin_id"]
   

    scar_fields = {}
    for scar, scar_sequence in check_scar.items():
        if scar_sequence:  # assuming non-empty string indicates presence
            name_entry = create_input(scar_frame, f"{scar}: Name (Required)")
            id_entry = create_input(scar_frame, f"{scar}: ID (Required)")
            description_entry = create_input(scar_frame, f"{scar}: Description (Optional)")
            scar_fields[scar] = (name_entry, id_entry, description_entry, scar_sequence)
    
    # Function to process data
    
    

    def submit_data():
        nonlocal data, scar_data
        data = {
            "plasmid_name": plasmid_name.get(),
            "plasmid_id": plasmid_id.get(),
            "plasmid_description": plasmid_description.get("1.0", tk.END).strip(),
            "cargo_name": cargo_name.get(),
            "cargo_id": cargo_id.get(),
            "cargo_description": cargo_description.get("1.0", tk.END).strip(),
            "marker_name": marker_name.get(),
            "marker_id": marker_id.get(),
            "marker_description": marker_description.get("1.0", tk.END).strip(),
            "origin_name": origin_name.get(),
            "origin_id": origin_id.get(),
            "origin_description": origin_description.get("1.0", tk.END).strip(),
        }

        for scar, (scar_name, scar_id, scar_description, scar_sequence) in scar_fields.items():
            name = scar_name.get()
            id = scar_id.get()
            description = scar_description.get()
            if name or id or description:  
                scar_data[scar] = [name, id, description, scar_sequence]
            if name == '' or id == '':
                messagebox.showerror("Error", "All requried fields must be filled")
                return
        print(scar_data)   
        if any(not data[field] for field in required_fields):
            messagebox.showerror("Error", "All required fields must be filled")
            return
          
        root.quit()

    submit_button = ttk.Button(main_frame, text="Submit", command=submit_data)
    submit_button.grid(columnspan=2, pady=10)
    root.mainloop()
    if data and all(data.get(field) for field in required_fields):  # Ensures processing only happens with valid input
            cargo = get_module(sequence, data["cargo_name"], data["cargo_id"], "cargo", data["cargo_description"])
            marker = get_module(sequence, data["marker_name"], data["marker_id"], "marker", data["marker_description"])
            origin = get_module(sequence, data["origin_name"], data["origin_id"], "origin", data["origin_description"])

            generated_plasmid = Plasmid(data["plasmid_name"],data["plasmid_id"],data["plasmid_description"])
            generated_plasmid.structure = [cargo,marker,origin]


            if scar_data:
                for key in scar_data:
                    parts = key.split()  
                    position = parts[1]  
                    module = parts[2]  
                    scar_data[key].extend([module, position])

            if scar_data:
                for scar, (scar_name, scar_id, scar_description, scar_sequence, module, position) in scar_data.items():
                    if scar:
                        scar_part = Part(scar_name, scar_id, scar_sequence, "scar", scar_description)
                        generated_plasmid.insert_scar(scar_part, module, position)








            return generated_plasmid
    else:
        print("Submission incomplete")
    

def serialise_genbank_dict(dict, folder_name):
        """Serialises in genbank format: Need input file path"""
        if not isinstance(folder_name,str):
             raise TypeError("Folder name must be of type: String")
        
        current_direct = os.getcwd()
        path = os.path.join(current_direct, folder_name)
        os.makedirs(path)

        for unique_id, plasmid in dict.items():
            sequence = plasmid.get_sequence()
        

            sequence_r = SeqRecord(
                Seq(sequence),
                id=plasmid.unique_id,
                name=plasmid.name,
                description=plasmid.description,
            )
            print(sequence_r)
            data_dict = (
                plasmid.gather_data()
            )  # Obtains data dictionary and uses this to add features of each part
            print(data_dict)
            conversion = {
                "re": "misc_feature",  # Dictionary to convert the role into the GenBank accepted role
                "cargo": "misc_feature",
                "abr": "misc_feature",
                "ori": "rep_origin",
            }

            for part, attribute in data_dict.items():  # Role Conversion
                if attribute[2] in conversion:
                    attribute[2] = conversion[attribute[2]]

            sequence_r.annotations["molecule_type"] = (
                "DNA"  # Annotate the sequence to show that its DNA
            )

            # Going through the data dictionary an using it to annotate each part of the sequence
            for part, attribute in data_dict.items():
                feature = SeqFeature(
                    FeatureLocation(attribute[0], attribute[4]),
                    type=attribute[2],
                    qualifiers={
                        "gene": [part],  # Gene name
                        "locus_tag": [attribute[1]],  # A unique identifier for the gene
                        "note": [attribute[3]],  # Any additional notes
                    },
                )
                sequence_r.features.append(feature)

            # Write the SeqRecord to a GenBank file
            with open(f"{path}\{unique_id}.gb", "w") as output_handle:
                SeqIO.write(sequence_r, output_handle, "genbank")

        print(f"GenBank files saved at {path}")



def process_plasmid_list(plasmid_list,module):
    construction_dictionary = {}
    
    for donor, recipient, desired_plasmid_id in plasmid_list:

        donor_name = donor.name
        recipient_name = recipient.name

        construction_dictionary[desired_plasmid_id] = []
        if module == "cargo":
            donor_insert = donor.structure[0]
            donor_insert_length = len(donor_insert.get_sequence()) - 4  # 4 due to sticky ends
            donor_linearised_plasmid_length = len(donor.get_sequence()) - donor_insert_length

            recipient_removed_sequence = recipient.structure[0]
            recipient_removed_sequence_length = len(recipient_removed_sequence.get_sequence()) - 4
            recipient_linearlised_plasmid_length = len(recipient.get_sequence()) - recipient_removed_sequence_length

        if module == "marker":
            donor_insert = donor.structure[2]
            donor_insert_length = len(donor_insert.get_sequence()) - 9  # 4 due to sticky ends
            donor_linearised_plasmid_length = len(donor.get_sequence()) - donor_insert_length

            recipient_removed_sequence = recipient.structure[2]
            recipient_removed_sequence_length = len(recipient_removed_sequence.get_sequence()) - 9
            recipient_linearlised_plasmid_length = len(recipient.get_sequence()) - recipient_removed_sequence_length

        if module == "origin":
            donor_insert = donor.structure[4]
            donor_insert_length = len(donor_insert.get_sequence()) - 4  # 4 due to sticky ends
            donor_linearised_plasmid_length = len(donor.get_sequence()) - donor_insert_length

            recipient_removed_sequence = recipient.structure[4]
            recipient_removed_sequence_length = len(recipient_removed_sequence.get_sequence()) - 4
            recipient_linearlised_plasmid_length = len(recipient.get_sequence()) - recipient_removed_sequence_length

        construction_dictionary[desired_plasmid_id].extend([
            donor_name,
            donor_insert_length,
            donor_linearised_plasmid_length,
            recipient_name,
            recipient_linearlised_plasmid_length,
            recipient_removed_sequence_length
        ])

        

    return construction_dictionary




def process_module_plasmid_list(plasmid_list, module):
    construction_dictionary = {}
    
    for recipient, insert, desired_plasmid_id in plasmid_list:

        recipient_name = recipient.name
        insert_name = insert.name
        construction_dictionary[desired_plasmid_id] = []

        if module == "cargo":
            insert = insert
            donor_insert_length = len(insert.get_sequence()) - 4  # 4 due to sticky ends

            recipient_removed_sequence = recipient.structure[0]
            recipient_removed_sequence_length = len(recipient_removed_sequence.get_sequence()) - 4
            recipient_linearlised_plasmid_length = len(recipient.get_sequence()) - recipient_removed_sequence_length

        if module == "marker":
            insert = insert
            donor_insert_length = len(insert.get_sequence()) - 9  # 4 due to sticky ends

            recipient_removed_sequence = recipient.structure[2]
            recipient_removed_sequence_length = len(recipient_removed_sequence.get_sequence()) - 9
            recipient_linearlised_plasmid_length = len(recipient.get_sequence()) - recipient_removed_sequence_length

        if module == "origin":
            insert = insert
            donor_insert_length = len(insert.get_sequence()) - 4  # 4 due to sticky ends

            recipient_removed_sequence = recipient.structure[4]
            recipient_removed_sequence_length = len(recipient_removed_sequence.get_sequence()) - 4
            recipient_linearlised_plasmid_length = len(recipient.get_sequence()) - recipient_removed_sequence_length

        construction_dictionary[desired_plasmid_id].extend([
            insert_name,
            donor_insert_length,
            recipient_name,
            recipient_linearlised_plasmid_length,
            recipient_removed_sequence_length
        ])

    return construction_dictionary



def import_single_parameters(csv_filename, input_csv):
    

    with open('Opentron_Param_Single.csv', mode='r', encoding='utf-8-sig') as csvfile:
        reader = csv.reader(csvfile)
        rows = list(reader)
       
        
        # Insert the plasmid and fragment data into the CSV
        plasmid_pairs = list(input_csv.values())
        for i, values in enumerate(plasmid_pairs):
            donor_details = [values[0],values[1],values[2]]
            recipient_details = [values[3],values[4],values[5]]

            col_number =  (i // 8) * 2 + 1
            
            counter = i % 8
            reset_number = 8*counter + 2
            print(i,donor_details)
            print(reset_number)


            rows[reset_number][col_number] = donor_details[0]
            rows[reset_number+1][col_number] = donor_details[1]
            rows[reset_number+2][col_number] = donor_details[2]

            rows[reset_number][col_number + 1] = recipient_details[0]
            rows[reset_number+1][col_number + 1] = recipient_details[1]
            rows[reset_number+2][col_number + 1] = recipient_details[2]

        current_direct = os.getcwd()
        path = os.path.join(current_direct, "Protocol CSVs", "Digestion")
        if not os.path.exists(path):
            os.makedirs(path)

        with open(f"{path}\{csv_filename}_experimental_data_single.csv", mode='w', newline='', encoding='utf-8') as file:
            writer = csv.writer(file)
            writer.writerows(rows)
        os.startfile(f"{path}\{csv_filename}_experimental_data_single.csv")




def import_module_single_parameters(csv_filename, input_csv):

    with open('Opentron_Param_Module_Single.csv', mode='r', encoding='utf-8-sig') as csvfile:
        reader = csv.reader(csvfile)
        rows = list(reader)
       
        
        # Insert the plasmid and fragment data into the CSV
        plasmid_pairs = list(input_csv.values())
        for i, values in enumerate(plasmid_pairs):

            donor_details = [values[0],values[1]]
            recipient_details = [values[2],values[3],values[4]]

            col_number =  (i // 8) * 2 + 1
            
            counter = i % 8
            reset_number = 8*counter + 2
            print(i,donor_details)
            print(reset_number)


            rows[reset_number][col_number] = donor_details[0]
            rows[reset_number+1][col_number] = donor_details[1]

            rows[reset_number][col_number + 1] = recipient_details[0]
            rows[reset_number+1][col_number + 1] = recipient_details[1]
            rows[reset_number+2][col_number + 1] = recipient_details[2]

        current_direct = os.getcwd()
        path = os.path.join(current_direct, "Protocol CSVs", "Digestion")
        if not os.path.exists(path):
            os.makedirs(path)

        with open(f"{path}\{csv_filename}_experimental_data_single.csv", mode='w', newline='', encoding='utf-8') as file:
            writer = csv.writer(file)
            writer.writerows(rows)
        os.startfile(f"{path}\{csv_filename}_experimental_data_single.csv")








def get_single_parameters():


        
    data = []
    with open("Opentron_Param_Single_Inputted.csv", mode='r', encoding='utf-8-sig') as csvfile:
        reader = csv.reader(csvfile)
        plasmid_count = 0
        for row in reader:
        # Check if it's a "Plasmid Name" row
            if "Plasmid Name" in row[0]:
                # Count all non-empty entries except for the first one (which is the text "Plasmid Name")
                for item in row[1:]:
                    if item:  # This checks if the cell is not empty
                        plasmid_count += 1
            # print(plasmid_count)

    full_cols, remaining_cells = divmod(plasmid_count, 16)
    print(remaining_cells)
    print(full_cols)
    # print(remaining_cells)

    with open("Opentron_Param_Single_Inputted.csv", mode='r', encoding='utf-8-sig') as csvfile:
        reader = csv.reader(csvfile)
        rows = list(reader)
    


    if full_cols > 0:
        for q in range(1,(full_cols*2)+1):
            for i in range(0,8):
                data.append([rows[(8*i)+1][q],rows[(8*i)+2][q],rows[(8*i)+5][q],rows[(8*i)+6][q],rows[(8*i)+7][q]])

        for q in range((full_cols*2) +1,(full_cols*2) +3):
            for i in range(0,int(remaining_cells/2)):
                data.append([rows[(8*i)+1][q],rows[(8*i)+2][q],rows[(8*i)+5][q],rows[(8*i)+6][q],rows[(8*i)+7][q]])
    if full_cols == 0:
        for q in range(1,3):
            for i in range(0,int(remaining_cells/2)):
                data.append([rows[(8*i)+1][q],rows[(8*i)+2][q],rows[(8*i)+5][q],rows[(8*i)+6][q],rows[(8*i)+7][q]])

    

    data_list = data
    print(data_list)   #Data list does now not include frag lengths






def import_multi_parameters(csv_filename, input_csv):
    

    with open('Opentron_Param_Multi.csv', mode='r', encoding='utf-8-sig') as csvfile:
        reader = csv.reader(csvfile)
        rows = list(reader)

    base_row = 6  # Example starting row

    for i, values in enumerate(input_csv.values()):
        donor_details = [values[0], values[1], values[2]]
        recipient_details = [values[3], values[4], values[5]]
        
        col_number = (i // 8) * 2 + 1
        block = i % 8
        row_offset = 5 * block + 2
        current_row = base_row + row_offset

        # Assign values safely
        rows[current_row][col_number] = donor_details[0]
        rows[current_row+1][col_number] = donor_details[1]
        rows[current_row+2][col_number] = donor_details[2]

        rows[current_row][col_number + 1] = recipient_details[0]
        rows[current_row+1][col_number + 1] = recipient_details[1]
        rows[current_row+2][col_number + 1] = recipient_details[2]

        current_direct = os.getcwd()
        path = os.path.join(current_direct, "Protocol CSVs","Digestion")
        if not os.path.exists(path):
            os.makedirs(path)

        with open(f"{path}\{csv_filename}_experimental_data_multi.csv", mode='w', newline='', encoding='utf-8') as file:
            writer = csv.writer(file)
            writer.writerows(rows)
        os.startfile(f"{path}\{csv_filename}_experimental_data_multi.csv")


def import_module_multi_parameters(csv_filename, input_csv):
    

    with open('Opentron_Param_Module_Multi.csv', mode='r', encoding='utf-8-sig') as csvfile:
        reader = csv.reader(csvfile)
        rows = list(reader)

    base_row = 6  

    for i, values in enumerate(input_csv.values()):
        donor_details = [values[0], values[1]]
        recipient_details = [values[2], values[3], values[4]]
        
        col_number = (i // 8) * 2 + 1
        block = i % 8
        row_offset = 5 * block + 2
        current_row = base_row + row_offset

     
        rows[current_row][col_number] = donor_details[0]
        rows[current_row+1][col_number] = donor_details[1]

        rows[current_row][col_number + 1] = recipient_details[0]
        rows[current_row+1][col_number + 1] = recipient_details[1]
        rows[current_row+2][col_number + 1] = recipient_details[2]

        current_direct = os.getcwd()
        path = os.path.join(current_direct, "Protocol CSVs","Digestion")
        if not os.path.exists(path):
            os.makedirs(path)

        with open(f"{path}\{csv_filename}_experimental_data_multi.csv", mode='w', newline='', encoding='utf-8') as file:
            writer = csv.writer(file)
            writer.writerows(rows)
        os.startfile(f"{path}\{csv_filename}_experimental_data_multi.csv")







def get_multi_parameters():


    data = []
    with open("Opentron_Param_Multi_Inputted.csv", mode='r', encoding='utf-8-sig') as csvfile:
        reader = csv.reader(csvfile)
        plasmid_count = 0
        for row in reader:
        # Check if it's a "Plasmid Name" row
            if "Plasmid Name" in row[0]:
                # Count all non-empty entries except for the first one (which is the text "Plasmid Name")
                for item in row[1:]:
                    if item:  # This checks if the cell is not empty
                        plasmid_count += 1
            # print(plasmid_count)

    full_cols, remaining_cells = divmod(plasmid_count, 16)
    print(remaining_cells)
    print(full_cols)
    # print(remaining_cells)

    with open("Opentron_Param_Multi_Inputted.csv", mode='r', encoding='utf-8-sig') as csvfile:
        reader = csv.reader(csvfile)
        rows = list(reader)
    
    data.append([rows[0][1],rows[1][1],rows[2][1],rows[3][1],rows[4][1]])

    if full_cols > 0:
        for q in range(1,(full_cols*2)+1):
            for i in range(0,8):
                data.append([rows[(5*i)+7][q],rows[(5*i)+8][q]])

        for q in range((full_cols*2) +1,(full_cols*2) +3):
            for i in range(0,int(remaining_cells/2)):
                data.append([rows[(5*i)+7][q],rows[(5*i)+8][q]])
    
    if full_cols == 0:
        for q in range(1,3):
            for i in range(0,int(remaining_cells/2)):
                data.append([rows[(5*i)+7][q],rows[(5*i)+8][q]])

    


    data_list = data
    print(data_list)   #Data list does now not include frag lengths



def calculate_volumes(insert_length, insert_concentration, recipient_length, recipient_concentration, total_dna=100, ratio=3):
    # Average molecular weight of a base pair in DNA
    mw_bp = 660

    # Calculate the molecular weights of the insert and recipient
    mw_insert = insert_length * mw_bp
    mw_recipient = recipient_length * mw_bp

    # Solve for moles of recipient
    total_mass_grams = total_dna * 1e-9  # Convert ng to grams
    # R * (ratio * mw_insert + mw_recipient) = total_mass_grams
    moles_recipient = total_mass_grams / (ratio * mw_insert + mw_recipient)
    
    # Calculate moles of insert
    moles_insert = ratio * moles_recipient

    # Convert moles back to mass in ng
    mass_insert_ng = moles_insert * mw_insert * 1e9
    mass_recipient_ng = moles_recipient * mw_recipient * 1e9
    
    # Calculate the volume of each DNA solution to use
    volume_insert = mass_insert_ng / insert_concentration
    volume_recipient = mass_recipient_ng / recipient_concentration
    
    return volume_insert, volume_recipient







def import_parameters_ligation(csv_filename, input_csv, protocol_type):
    
    current_direct = os.getcwd()
    file = f"{csv_filename}_input_concentation_values.csv"
    path = os.path.join(current_direct, "Protocol CSVs", "Ligation",file)
    if os.path.exists(path):
        return None

    with open('Opentron_Param_Ligation.csv', mode='r', encoding='utf-8-sig') as csvfile:
        reader = csv.reader(csvfile)
        rows = list(reader)
       
        
        # Insert the plasmid and fragment data into the CSV
        plasmid_pairs = list(input_csv.values())
        for i, values in enumerate(plasmid_pairs):
            if protocol_type == "plasmid":
                donor_details = [values[0],values[1]]
                recipient_details = [values[3],values[4]]
            elif protocol_type == "module":
                donor_details = [values[0],values[1]]
                recipient_details = [values[2],values[4]]

            col_number =  (i // 8) * 2 + 1
            
            counter = i % 8
            reset_number = 6*counter + 2
            print(i,donor_details)
            print(reset_number)


            rows[reset_number][col_number] = donor_details[0]
            rows[reset_number+1][col_number] = donor_details[1]

            rows[reset_number][col_number + 1] = recipient_details[0]
            rows[reset_number+1][col_number + 1] = recipient_details[1]

        current_direct = os.getcwd()
        path = os.path.join(current_direct, "Protocol CSVs", "Ligation")
        if not os.path.exists(path):
            os.makedirs(path)

        with open(f"{path}\{csv_filename}_input_concentation_values.csv", mode='w', newline='', encoding='utf-8') as file:
            writer = csv.writer(file)
            writer.writerows(rows)


def generate_ligation_experiment(filename):
    
    current_direct = os.getcwd()
    path = os.path.join(current_direct, "Protocol CSVs", "Ligation")
        
    data = []
    with open(f"{path}\{filename}_input_concentation_values.csv", mode='r', encoding='utf-8-sig') as csvfile:
        reader = csv.reader(csvfile)
        plasmid_count = 0
        for row in reader:
        # Check if it's a "Plasmid Name" row
            if "Plasmid Name" in row[0]:
                # Count all non-empty entries except for the first one (which is the text "Plasmid Name")
                for item in row[1:]:
                    if item:  # This checks if the cell is not empty
                        plasmid_count += 1
            # print(plasmid_count)

    full_cols, remaining_cells = divmod(plasmid_count, 16)

    with open(f"{path}\{filename}_input_concentation_values.csv", mode='r', encoding='utf-8-sig') as csvfile:
        reader = csv.reader(csvfile)
        rows = list(reader)
    


    if full_cols > 0:
        for q in range(1,(full_cols*2)+1):
            for i in range(0,8):
                data.append([rows[(6*i)+1][q],rows[(6*i)+2][q],rows[(6*i)+3][q],rows[(6*i)+4][q]])

        for q in range((full_cols*2) +1,(full_cols*2) +3):
            for i in range(0,int(remaining_cells/2)):
                data.append([rows[(6*i)+1][q],rows[(6*i)+2][q],rows[(6*i)+3][q],rows[(6*i)+4][q]])
    if full_cols == 0:
        for q in range(1,3):
            for i in range(0,int(remaining_cells/2)):
                data.append([rows[(6*i)+1][q],rows[(6*i)+2][q],rows[(6*i)+3][q],rows[(6*i)+4][q]])

    

    data_list = data
    
    def sort_custom_simple(data):
    # Simple sorting function
        sorted_data = sorted(data, key=lambda x: (x[0][0], int(x[0][1:])))
        return sorted_data
 
    sorted_data = sort_custom_simple(data_list)

    def group_in_pairs(data):
        return [data[i:i+2] for i in range(0, len(data), 2)]

    paired_data = group_in_pairs(sorted_data)
    for i in paired_data:
        # if i[0][2] or i[0][3] or i[1][2] or i[1][3] == '':
        #     raise ValueError("Please make sure all concentrations are filled in")
        
        vol_insert, vol_recipient = calculate_volumes(float(i[0][2]), float(i[0][3]), float(i[1][2]), float(i[1][3]))
        if vol_insert + vol_recipient <=18:
            i[0].append(round(vol_insert,2))
            i[1].append(round(vol_recipient,2))
        elif vol_insert + vol_recipient > 18:
            i[0].append(0)
            i[1].append(0)

   
    with open(f"{path}\{filename}_input_concentation_values.csv", mode='r', encoding='utf-8-sig') as csvfile:
        reader = csv.reader(csvfile)
        data = list(reader)

    # Transpose the data (convert rows to columns)
    # Zip the rows to create columns, then convert each to a list
        transposed_data = [list(col) for col in zip(*data)]
     
    def find_index_in_nested_list(data, target):
        for row_index, row in enumerate(data):
            if target in row:
                return row_index , row.index(target) 
        return None  # Return None if the target is not found
    
    for i in paired_data:

        col1, row1 = find_index_in_nested_list(transposed_data,i[0][0])
        transposed_data[col1][row1+4] = i[0][4]

        col2, row2 = find_index_in_nested_list(transposed_data,i[1][0])
        transposed_data[col2][row2+4] = i[1][4]
 
    np_data = np.array(transposed_data)

    row_data = np_data.T

    row_list = row_data.tolist()



    with open(f"{path}\{filename}_experimental_data.csv", mode='w', newline='', encoding='utf-8') as file:
            writer = csv.writer(file)
            writer.writerows(row_list)
    os.startfile(f"{path}\{filename}_experimental_data.csv")



  
def import_insert_parameters(csv_filename, input_csv):

    with open('Opentron_Param_Insert.csv', mode='r', encoding='utf-8-sig') as csvfile:
        reader = csv.reader(csvfile)
        rows = list(reader)
       
        
        # Insert the plasmid and fragment data into the CSV
        # plasmid_pairs = list(input_csv.values())
        for i, values in enumerate(input_csv):

            donor_details = [values[0],values[1]]
            recipient_details = [values[2],values[3],values[4]]
            enzyme_details = [values[5],values[6]]

            col_number =  (i // 8) * 2 + 1
            
            counter = i % 8
            reset_number = 10*counter + 2
            print(i,donor_details)
            print(reset_number)


            rows[reset_number][col_number] = donor_details[0]
            rows[reset_number+1][col_number] = donor_details[1]
            rows[reset_number+6][col_number] = enzyme_details[0]
            rows[reset_number+7][col_number] = enzyme_details[1]


            rows[reset_number][col_number + 1] = recipient_details[0]
            rows[reset_number+1][col_number + 1] = recipient_details[1]
            rows[reset_number+2][col_number + 1] = recipient_details[2]
            rows[reset_number+6][col_number+1] = enzyme_details[0]
            rows[reset_number+7][col_number+1] = enzyme_details[1]

        current_direct = os.getcwd()
        path = os.path.join(current_direct, "Protocol CSVs", "Digestion")
        if not os.path.exists(path):
            os.makedirs(path)

        with open(f"{path}\{csv_filename}_experimental_data.csv", mode='w', newline='', encoding='utf-8') as file:
            writer = csv.writer(file)
            writer.writerows(rows)
        os.startfile(f"{path}\{csv_filename}_experimental_data.csv")


input1 = [
    ['tester1', 289, 'plasmid1test', 16, 2086, 'BamHI', 'SalI'],
    ['tester1', 289, 'plasmid2test', 16, 2086, 'EcoRI', 'PstI'],
    ['tester2', 345, 'plasmid3test', 150, 1999, 'SacI', 'HindIII'],
    ['tester3', 410, 'plasmid4test', 200, 1800, 'BamHI', 'SalI'],
    ['tester4', 512, 'plasmid5test', 300, 1700, 'EcoRI', 'HindIII'],
    ['tester5', 620, 'plasmid6test', 400, 1600, 'SacI', 'PstI'],
    ['tester6', 730, 'plasmid7test', 500, 1500, 'BamHI', 'SalI'],
    ['tester7', 840, 'plasmid8test', 600, 1400, 'EcoRI', 'PstI'],
    ['tester8', 950, 'plasmid9test', 700, 1300, 'SacI', 'HindIII'],
    ['tester9', 1060, 'plasmid10test', 800, 1200, 'BamHI', 'SalI'],
    ['tester10', 1170, 'plasmid11test', 900, 1100, 'EcoRI', 'HindIII'],
    ['tester11', 1280, 'plasmid12test', 1000, 1000, 'SacI', 'PstI'],
    ['tester12', 1390, 'plasmid13test', 1100, 900, 'BamHI', 'SalI'],
    ['tester13', 1500, 'plasmid14test', 1200, 800, 'EcoRI', 'HindIII'],
    ['tester14', 1610, 'plasmid15test', 1300, 700, 'SacI', 'PstI'],
    ['tester15', 1720, 'plasmid16test', 1400, 600, 'BamHI', 'SalI'],
    ['tester16', 1830, 'plasmid17test', 1500, 500, 'EcoRI', 'HindIII'],
    ['tester17', 1940, 'plasmid18test', 1600, 400, 'SacI', 'PstI'],
    ['tester18', 2050, 'plasmid19test', 1700, 300, 'BamHI', 'SalI'],
    ['tester19', 2160, 'plasmid20test', 1800, 200, 'EcoRI', 'HindIII']
]
