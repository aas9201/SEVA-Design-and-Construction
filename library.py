import sqlite3
import json
from plasmid import Part,Module,Plasmid
import re
from utils import sort_sites, get_module, get_scars
import os
import tkinter as tk
from tkinter import ttk, messagebox
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio import SeqIO
import tkinter as tk
from tkinter import ttk, messagebox
from sbol_serialisation import serialise_sbol

class Library:
    
    def __init__(self, database_name):
        self.database_name = database_name
        self.db_path = self._get_full_db_path()
        self._initialise_database()
    

    def _get_full_db_path(self):
        # Ensure the directory exists
        current_dict = os.getcwd()
        database_folder = os.path.join(current_dict,"Database Folder")
        if not os.path.exists(database_folder):
            os.makedirs(database_folder)
        # Combine the directory path with the database file name
        return os.path.join(database_folder, self.database_name + ".db")
    

    def _initialise_database(self):
        with sqlite3.connect(self.db_path) as conn:
            cursor = conn.cursor()
            cursor.execute('''
                CREATE TABLE IF NOT EXISTS plasmids (
                    unique_id TEXT PRIMARY KEY,
                    plasmid_data TEXT NOT NULL
                )
            ''')
            cursor.execute('''
                CREATE TABLE IF NOT EXISTS modules (
                    unique_id TEXT PRIMARY KEY,
                    module_data TEXT NOT NULL
                )
            ''')
            cursor.execute('''
                CREATE TABLE IF NOT EXISTS parts (
                    unique_id TEXT PRIMARY KEY,
                    part_data TEXT NOT NULL
                )
            ''')
            cursor.execute('''
                CREATE TABLE IF NOT EXISTS inserts (
                    unique_id TEXT PRIMARY KEY,
                    insert_data TEXT NOT NULL
                )
            ''')
            conn.commit()

    def _save_plasmid(self,plasmid):
        if not isinstance(plasmid, Plasmid):
            raise TypeError("Only plasmid objects can be added")
        plasmid_data = json.dumps(plasmid.to_dict(), default=str)
        try:
            with sqlite3.connect(self.db_path) as conn:
                cursor = conn.cursor()
                cursor.execute('INSERT INTO plasmids (unique_id, plasmid_data) VALUES (?, ?)', 
                            (plasmid.unique_id, plasmid_data))
                conn.commit()
        except sqlite3.DatabaseError as e:
            print(f"Database error: {e}")

    def _save_insert(self,insert):
        insert_data = json.dumps(insert, default=str)
        try:
            with sqlite3.connect(self.db_path) as conn:
                cursor = conn.cursor()
                cursor.execute('INSERT INTO inserts (unique_id, insert_data) VALUES (?, ?)', 
                            (insert["unique_id"], insert_data))
                conn.commit()
        except sqlite3.DatabaseError as e:
            print(f"Database error: {e}")


    def _save_module(self, module):
        if not isinstance(module, Module):
            raise TypeError("Only Module objects can be added")
        module_data = json.dumps(module.to_dict(), default=str)
        with sqlite3.connect(self.db_path) as conn:
            cursor = conn.cursor()
            cursor.execute('INSERT INTO modules (unique_id, module_data) VALUES (?, ?)', 
                        (module.unique_id, module_data))
            conn.commit()

    def _save_part(self, part):
        if not isinstance(part, Part):
            raise TypeError("Only Part objects can be added")
        part_data = json.dumps(part.to_dict(), default=str)
        with sqlite3.connect(self.db_path) as conn:
            cursor = conn.cursor()
            cursor.execute('INSERT INTO parts (unique_id, part_data) VALUES (?, ?)', 
                        (part.unique_id, part_data))
            conn.commit()



    def _plasmid_exists(self, unique_id):
        with sqlite3.connect(self.db_path) as conn:
            cursor = conn.cursor()
            cursor.execute('SELECT plasmid_data FROM plasmids WHERE unique_id = ?', (unique_id,))
            plasmid_data = cursor.fetchone()

        if plasmid_data:
            return True
        else:
            return False
        
    def _insert_exists(self, unique_id):
        with sqlite3.connect(self.db_path) as conn:
            cursor = conn.cursor()
            cursor.execute('SELECT insert_data FROM inserts WHERE unique_id = ?', (unique_id,))
            insert_data = cursor.fetchone()

        if insert_data:
            return True
        else:
            return False
        
    def _module_exists(self, unique_id):
        with sqlite3.connect(self.db_path) as conn:
            cursor = conn.cursor()
            cursor.execute('SELECT module_data FROM modules WHERE unique_id = ?', (unique_id,))
            module_data = cursor.fetchone()

        if module_data:
            return True
        else:
            return False
        
    def _part_exists(self, unique_id):
        with sqlite3.connect(self.db_path) as conn:
            cursor = conn.cursor()
            cursor.execute('SELECT part_data FROM parts WHERE unique_id = ?', (unique_id,))
            part_data = cursor.fetchone()

        if part_data:
            return True
        else:
            return False


    def add_plasmid(self, *plasmids):
        for plasmid in plasmids:
            existing_plasmid = self._plasmid_exists(plasmid.unique_id)
            if existing_plasmid:
                print(f"Plasmid with Unique ID: {plasmid.unique_id} is already within database.")
                continue  # Skip to next plasmidd  
            else:
                for module in plasmid.structure:
                    if module._dynamic:  # Only add dynamic modules
                        for part in module.structure:
                            if part._dynamic:  # Only add dynamic parts
                                self.add_part(part)
                        self.add_module(module)
            
            self._save_plasmid(plasmid)
            print(f"Plasmid: '{plasmid.unique_id}', saved and processed.")

    


    def _constructor_add_plasmid(self, plasmid):
        existing_plasmid = self._plasmid_exists(plasmid.unique_id)
        if existing_plasmid:
            print(f"Plasmid with Unique ID: {plasmid.unique_id} is already within database.")
            return False  
        else:
            for module in plasmid.structure:
                for part in module.structure:
                    self.add_part(part)
                self.add_module(module)

            self._save_plasmid(plasmid)
            print(f"Plasmid: '{plasmid.unique_id}', saved and processed.")
            return True

    def _constructor_add_insert(self, insert):
        existing_plasmid = self._insert_exists(insert["unique_id"])
        if existing_plasmid:
            print(f"Insert with Unique ID: {insert['unique_id']} is already within database.")
            return False  
        else:
            self._save_insert(insert)
            print(f"Insert: {insert['unique_id']}, saved and processed.")
            return True

    def add_module(self, module):
        if not module._dynamic:
            return  # Skip saving non-dynamic modules
        existing_module = self._module_exists(module.unique_id)
        if existing_module:
            print(f"Module with Unique ID: {module.unique_id} is already within database.")
            return existing_module 
        else:
            self._save_module(module)
            print(f"Module: '{module.unique_id}', saved and processed.")


    def _constructor_add_module(self, module):
        existing_module = self._module_exists(module.unique_id)
        if existing_module:
            print(f"Plasmid with Unique ID: {module.unique_id} is already within database.")
            return False  
        
        else:
            self._save_module(module)
            print(f"Module: '{module.unique_id}', saved and processed.")
            return True
        

    def add_part(self, part):
        existing_part = self._part_exists(part.unique_id)
        if existing_part:
            print(f"Part with Unique ID: {part.unique_id} is already within database.")
            return existing_part  
        else:
            self._save_part(part)
            print(f"Part: '{part.unique_id}', saved and processed.")



    def add_plasmid_dict(self,plasmid_dict):

        for unique_id, plasmid_object in plasmid_dict.items():
            self.add_plasmid(plasmid_object)


    def delete_plasmid(self, unique_id):
        try:
            with sqlite3.connect(self.db_path) as conn:
                cursor = conn.cursor()
            
                cursor.execute('DELETE FROM plasmids WHERE unique_id = ?', (unique_id,))
                conn.commit() 
            if cursor.rowcount == 0: 
                print(f"No plasmid found with the ID: {unique_id}")
            else:
                print(f"Plasmid with unique_id {unique_id} deleted successfully.")
        except sqlite3.DatabaseError as e:
            print(f"Database error: {e}")


    def delete_module(self, unique_id):
        try:
            with sqlite3.connect(self.db_path) as conn:
                cursor = conn.cursor()
               
                cursor.execute('DELETE FROM modules WHERE unique_id = ?', (unique_id,))
                conn.commit()  
            if cursor.rowcount == 0:  
                print(f"No module found with the ID: {unique_id}")
            else:
                print(f"Module with unique_id {unique_id} deleted successfully.")
        except sqlite3.DatabaseError as e:
            print(f"Database error: {e}")


    def delete_part(self, unique_id):
        try:
            with sqlite3.connect(self.db_path) as conn:
                cursor = conn.cursor()
               
                cursor.execute('DELETE FROM parts WHERE unique_id = ?', (unique_id,))
                conn.commit() 
            if cursor.rowcount == 0:  
                print(f"No part found with the ID: {unique_id}")
            else:
                print(f"Part with unique_id {unique_id} deleted successfully.")
        except sqlite3.DatabaseError as e:
            print(f"Database error: {e}")



    def print_plasmid_table(self):
        try:
            with sqlite3.connect(self.db_path) as conn:
                cursor = conn.cursor()
                cursor.execute('SELECT unique_id, plasmid_data FROM plasmids')  
                rows = cursor.fetchall()  
            if rows:
                print("Unique ID | Plasmid Data")
                for row in rows:
                    print(f"{row[0]} | {row[1]}")
            else:
                print("No data found in the plasmids table.")
        except sqlite3.DatabaseError as e:
            print(f"Database error: {e}")


    def print_modules_table(self):
        try:
            with sqlite3.connect(self.db_path) as conn:
                cursor = conn.cursor()
                cursor.execute('SELECT unique_id, module_data FROM modules') 
                rows = cursor.fetchall() 
            if rows:
                print("Unique ID | Module Data")
                for row in rows:
                    print(f"{row[0]} | {row[1]}")
            else:
                print("No data found in the modules table.")
        except sqlite3.DatabaseError as e:
            print(f"Database error: {e}")


    def print_parts_table(self):
        try:
            with sqlite3.connect(self.db_path) as conn:
                cursor = conn.cursor()
                cursor.execute('SELECT unique_id, part_data FROM parts')  
                rows = cursor.fetchall()  
            if rows:
                print("Unique ID | Part Data")
                for row in rows:
                    print(f"{row[0]} | {row[1]}")
            else:
                print("No data found in the parts table.")
        except sqlite3.DatabaseError as e:
            print(f"Database error: {e}")


    def _retrieve_all_json_plasmids(self):
        with sqlite3.connect(self.db_path) as conn:
            cursor = conn.cursor()
            cursor.execute('SELECT unique_id, plasmid_data FROM plasmids')  
            data = cursor.fetchall() 
        return data
    

    def _retrieve_all_json_modules(self):
        with sqlite3.connect(self.db_path) as conn:
            cursor = conn.cursor()
            cursor.execute('SELECT unique_id, module_data FROM modules') 
            data = cursor.fetchall()  
        return data
    
    def _retrieve_all_json_inserts(self):
        with sqlite3.connect(self.db_path) as conn:
            cursor = conn.cursor()
            cursor.execute('SELECT unique_id, insert_data FROM inserts') 
            data = cursor.fetchall()  
        return data
    

    def _retrieve_all_json_parts(self):
        with sqlite3.connect(self.db_path) as conn:
            cursor = conn.cursor()
            cursor.execute('SELECT unique_id, part_data FROM parts')  
            data = cursor.fetchall()  
        return data


    def _retrieve_one_json_plasmid(self, unique_id):
        with sqlite3.connect(self.db_path) as conn:
            cursor = conn.cursor()
            cursor.execute('SELECT plasmid_data FROM plasmids WHERE unique_id = ?', (unique_id,))
            plasmid_data = cursor.fetchone()

        if plasmid_data:
            return plasmid_data[0]
        else:
            print(f"Plasmid with Unique ID '{unique_id}' not found")
    

    def _retrieve_one_json_module(self, unique_id):
        with sqlite3.connect(self.db_path) as conn:
            cursor = conn.cursor()
            cursor.execute('SELECT module_data FROM modules WHERE unique_id = ?', (unique_id,))
            module_data = cursor.fetchone()

        if module_data:
            return module_data[0]
        else:
            print(f"Module with Unique ID '{unique_id}' not found")


    def _retrieve_one_json_part(self, unique_id):
        with sqlite3.connect(self.db_path) as conn:
            cursor = conn.cursor()
            cursor.execute('SELECT part_data FROM parts WHERE unique_id = ?', (unique_id,))
            part_data = cursor.fetchone()

        if part_data:
            return part_data[0]
        else:
            print(f"Part with Unique ID '{unique_id}' not found")

    def retrieve_insert(self, unique_id):
        with sqlite3.connect(self.db_path) as conn:
            cursor = conn.cursor()
            cursor.execute('SELECT insert_data FROM inserts WHERE unique_id = ?', (unique_id,))
            insert_json_data = cursor.fetchone()

        if insert_json_data:
            # Access the first element of the tuple which contains the JSON data
            insert_data = json.loads(insert_json_data[0])  # Assuming the JSON data is in the first column
            return insert_data
        else:
            print(f"Part with Unique ID '{unique_id}' not found")
            return None



    def retrieve_plasmid(self, *unique_ids):
        plasmids = []
        for unique_id in unique_ids:
            json_plasmid_data = self._retrieve_one_json_plasmid(unique_id)
            plasmid_data = json.loads(json_plasmid_data)
            plasmid = self.json_to_plasmid(plasmid_data)
            plasmids.append(plasmid)
        
        if len(plasmids) == 1:
            return plasmids[0]  
        return plasmids  
    
    def retrieve_module(self, *unique_ids):
        modules = []
        for unique_id in unique_ids:
            json_module_data = self._retrieve_one_json_module(unique_id)
            module_data = json.loads(json_module_data)
            module = self.json_to_module(module_data)
            modules.append(module)
        
        if len(modules) == 1:
            return modules[0] 
        return modules 

    def retrieve_part(self, *unique_ids):
        parts = []
        for unique_id in unique_ids:
            json_part_data = self._retrieve_one_json_part(unique_id)
            part_data = json.loads(json_part_data)
            part = self.json_to_part(part_data)
            parts.append(part)
        
        if len(parts) == 1:
            return parts[0]  
        return parts  

    def retrieve_part(self, *unique_ids):
        parts = []
        for unique_id in unique_ids:
            json_part_data = self._retrieve_one_json_part(unique_id)
            part_data = json.loads(json_part_data)
            part = self.json_to_part(part_data)
            parts.append(part)
        
        if len(parts) == 1:
            return parts[0]  
        return parts  


    def json_to_plasmid(self,plasmid_data):
        modules = []
        for module_data in plasmid_data['structure']:  
            parts = []
            for part_data in module_data['structure']:  
                part = Part(
                    name=part_data['name'],
                    unique_id=part_data['unique_id'],
                    sequence=part_data['sequence'],
                    role=part_data['role'],
                    description=part_data['description']
                )
                parts.append(part)
            module = Module(
                name=module_data['name'],
                unique_id=module_data['unique_id'],
                module_type=module_data['module_type'],
                description=module_data['description'],
                structure=parts
            )
            modules.append(module)

     
        plasmid = Plasmid(
            name=plasmid_data['name'],
            unique_id=str(plasmid_data['unique_id']),
            description=plasmid_data['description']
        )
        plasmid.structure = modules

        if 'scars' in plasmid_data:
            for scar_id, scar_info in plasmid_data['scars'].items():
                scar_part = Part(
                    name=scar_info['name'],
                    unique_id=scar_id,
                    sequence=scar_info['sequence'],
                    role='scar',
                    description=scar_info['description']
                )
                plasmid.insert_scar(scar_part, scar_info['module'], scar_info['position'])

        return plasmid 



    def json_to_module(self,module_data):
        parts = []
        for part_data in module_data['structure']:  
            part = Part(
            name=part_data['name'],
            unique_id=part_data['unique_id'],
            sequence=part_data['sequence'],
            role=part_data['role'],
            description=part_data['description']
            )
            parts.append(part)
        module = Module(
        name=module_data['name'],
        unique_id=module_data['unique_id'],
        module_type=module_data['module_type'],
        description=module_data['description'],
        structure=parts
        )
        return module

    def json_to_part(self,part_data):
        part = Part(
    name=part_data['name'],
    unique_id=part_data['unique_id'],
    sequence=part_data['sequence'],
    role=part_data['role'],
    description=part_data['description']
    )
        return part

    def retrieve_all_plasmids(self):
        json_data  = self._retrieve_all_json_plasmids()
        new_dict = {}
        for unique_id, json_convert in json_data:
            plasmid_data = json.loads(json_convert)
            plasmid = self.json_to_plasmid(plasmid_data)
            
            new_dict[unique_id] = plasmid

        return new_dict
    

    def retrieve_all_modules(self):
        json_data  = self._retrieve_all_json_modules()
        new_dict = {}
        for unique_id, json_convert in json_data:
            module_data = json.loads(json_convert)
            module = self.json_to_module(module_data)
            
            new_dict[unique_id] = module

        return new_dict
    

    def retrieve_all_parts(self):
        json_data  = self._retrieve_all_json_parts()
        new_dict = {}
        for unique_id, json_convert in json_data:
            part_data = json.loads(json_convert)
            part = self.json_to_part(part_data)
            
            new_dict[unique_id] = part

        return new_dict


    def save_json_file(self, file_name):
        current_direct = os.getcwd()
        path = os.path.join(current_direct, "Plasmid Storage")
        if not os.path.exists(path):
            os.makedirs(path)

        file_path = os.path.join(path, file_name)
        plasmid_data = self._retrieve_all_json_plasmids()

        with open(file_path, 'w') as file:
            json.dump(plasmid_data, file, indent=3)


    def load_file(self, file_name):
        current_direct = os.getcwd()
        path = os.path.join(current_direct, "Plasmid Storage")
        file_path = os.path.join(path, file_name)

        new_dict = {}
        with open(file_path, 'r') as file:
            plasmid_data = json.load(file)
        for unique_id, json_convert in plasmid_data:
            plasmid_data = json.loads(json_convert)
            plasmid = self.json_to_plasmid(plasmid_data)
            
            new_dict[unique_id] = plasmid

        return new_dict
    
    def _retrieve_all_module_id(self):
        data = self._retrieve_all_json_plasmids()
        module_id_full = []  
        for unique_id, plasmid_data in data:
            plasmids = json.loads(plasmid_data)
            module_ids = []  
            for module_data in plasmids['structure']:
                module_id = module_data["unique_id"]
                module_ids.append(module_id)
            module_id_full.append(module_ids)  

        return module_id_full

    def _module_to_plasmid(self):
        #makes dictionary with keys being module ids in a plasmid with value as plasmid ID
        module_to_plasmid = {}
        data = self._retrieve_all_json_plasmids()

        for plasmid_id, plasmid_data in data:
            plasmids = json.loads(plasmid_data)
            module_ids = tuple((module["unique_id"] for module in plasmids["structure"]))
            module_to_plasmid[module_ids] = plasmid_id
        
        return module_to_plasmid
    


    def serialise_genbank(self, folder_name, *plasmid_ids):
        """Serialises in genbank format: Need input file path"""
        if not isinstance(folder_name,str):
             raise TypeError("Folder name must be of type: String")
        
        current_direct = os.getcwd()
        path = os.path.join(current_direct, folder_name)
        if not os.path.exists(path):
            os.makedirs(path)

        plasmid_list = self.retrieve_plasmid(*plasmid_ids)

        if not isinstance(plasmid_list, list):
            plasmid_list = [plasmid_list]
        

        for plasmid in plasmid_list:
            sequence = plasmid.get_sequence()

        

            sequence_r = SeqRecord(
                Seq(sequence),
                id=plasmid.unique_id,
                name=plasmid.name,
                description=plasmid.description,
            )

            data_dict = (
                plasmid.gather_data()
            )  # Obtains data dictionary and uses this to add features of each part
    
            conversion = {
                "re": "misc_feature",  # Dictionary to convert the role into the GenBank accepted role
                "cargo": "misc_feature",
                "abr": "misc_feature",
                "ori": "rep_origin",
                "orit": "oriT",
                "gene": "gene",
                "enhancer": "enhancer",
                "attenuator": "attenuator",
                "misc":"misc_feature",
                "scar": "misc_feature",
                "operator": "regulatory",
                "rbs":"RBS"
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
            with open(f"{path}\{plasmid.unique_id}.gb", "w") as output_handle:
                SeqIO.write(sequence_r, output_handle, "genbank")

        print(f"GenBank files saved at {path}")

    def serialise_fasta(self, folder_name, *plasmid_ids):
        """Serialises in FASTA format: Need input file path"""
        if not isinstance(folder_name,str):
             raise TypeError("Folder name must be of type: String")
        
        current_direct = os.getcwd()
        path = os.path.join(current_direct, folder_name)
        if not os.path.exists(path):
            os.makedirs(path)

        plasmid_list = self.retrieve_plasmid(*plasmid_ids)

        if not isinstance(plasmid_list, list):
            plasmid_list = [plasmid_list]
        

        for plasmid in plasmid_list:
            sequence = plasmid.get_sequence()

            sequence_r = SeqRecord(
                Seq(sequence),
                id=plasmid.unique_id,
                name=plasmid.name,
                description=plasmid.description,
            )
            with open(f"{path}\{plasmid.unique_id}.fasta", "w") as output_handle:
                SeqIO.write(sequence_r, output_handle, "fasta")

    def import_plasmid(self, sequence):

       
        check_scar = get_scars(sequence)
        data = None
        scar_data = {}
        self.root = tk.Tk()
        self.root.style = ttk.Style()
        self.root.style.theme_use('clam')
        self.root.style.theme_use('clam')  # Modern theme
        self.root.style.configure('TFrame', background='#f0f0f0')
        self.root.style.configure('TButton', font=('Helvetica', 10))
        self.root.style.configure('TLabel', background='#f0f0f0', font=('Helvetica', 10))
        self.root.style.configure('Treeview', background='white', foreground='black', fieldbackground='white')
        self.root.style.map('Treeview', background=[('selected', '#0078d7')])
        self.root.title("Plasmid Module Details")

       
        main_frame = ttk.Frame(self.root, padding="10")
        main_frame.pack(fill='both', expand=True)

        plasmid_frame = ttk.Frame(main_frame, padding="10")
        plasmid_frame.grid(row=0, column=0, sticky='nswe')

        
        scar_frame = ttk.Frame(main_frame, padding="10")
        scar_frame.grid(row=0, column=1, sticky='nswe', padx=(10, 0))



      
        def create_input(container, label, is_text=False):
            ttk.Label(container, text=label).grid(columnspan=2, sticky='w')
            if is_text:
                input_widget = tk.Text(container, height=3, width=40, font = ("Helvetica", 12))
            else:
                input_widget = ttk.Entry(container, width=40, font = ("Helvetica", 12))
            input_widget.grid(columnspan=2, sticky='ew')
            return input_widget

       
        plasmid_name = create_input(plasmid_frame, "Plasmid Name (Required)")
        plasmid_id = create_input(plasmid_frame, "Plasmid ID (Required)")
        plasmid_description = create_input(plasmid_frame, "Plasmid Description (Not Required)", is_text=True)

        cargo_name = create_input(plasmid_frame, "Cargo Name (Required)")
        cargo_id = create_input(plasmid_frame, "Cargo ID (Required)")
        cargo_description = create_input(plasmid_frame, "Cargo Description (Not Required)", is_text=True)

        marker_name = create_input(plasmid_frame, "Marker Name (Required)")
        marker_id = create_input(plasmid_frame, "Marker ID (Required)")
        marker_description = create_input(plasmid_frame, "Marker Description (Not Required)", is_text=True)

        replication_name = create_input(plasmid_frame, "Replication Name (Required)")
        replication_id = create_input(plasmid_frame, "Replication ID (Required)")
        replication_description = create_input(plasmid_frame, "Replication Description (Not Required)", is_text=True)
        required_fields = ["plasmid_name", "plasmid_id", "cargo_name", "cargo_id", "marker_name", "marker_id", "replication_name", "replication_id"]
    

        scar_fields = {}
        for scar, scar_sequence in check_scar.items():
            if scar_sequence: 
                name_entry = create_input(scar_frame, f"{scar}: Name (Required)")
                id_entry = create_input(scar_frame, f"{scar}: ID (Required)")
                description_entry = create_input(scar_frame, f"{scar}: Description (Optional)")
                scar_fields[scar] = (name_entry, id_entry, description_entry, scar_sequence)
        
        
        

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
                "replication_name": replication_name.get(),
                "replication_id": replication_id.get(),
                "replication_description": replication_description.get("1.0", tk.END).strip(),
            }
            print(scar_fields)
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
            
            self.root.quit()

        submit_button = ttk.Button(main_frame, text="Submit", command=submit_data)
        submit_button.grid(columnspan=2, pady=10)
        self.root.mainloop()
        if data and all(data.get(field) for field in required_fields):  # Ensures processing only happens with valid input
                cargo = get_module(sequence, data["cargo_name"], data["cargo_id"], "cargo", data["cargo_description"])
                marker = get_module(sequence, data["marker_name"], data["marker_id"], "marker", data["marker_description"])
                replication = get_module(sequence, data["replication_name"], data["replication_id"], "replication", data["replication_description"])

                generated_plasmid = Plasmid(data["plasmid_name"],data["plasmid_id"],data["plasmid_description"])
                generated_plasmid.structure = [cargo,marker,replication]


                if scar_data:
                    for key in scar_data:
                        parts = key.split()  
                        position = parts[1]  
                        module = parts[2]  
                        scar_data[key].extend([module, position])

                    for scar, (scar_name, scar_id, scar_description, scar_sequence, module, position) in scar_data.items():
                        if scar:
                            scar_part = Part(scar_name, scar_id, scar_sequence, "scar", scar_description)
                            generated_plasmid.insert_scar(scar_part, module, position)


                reference_point = sequence.find(
                "TTAATTAA"
            )  # This just makes PacI as the starting reference point of every plasmid sequence due to it being circular
       
                ordered_sequence = (
                    sequence[reference_point:] + sequence[:reference_point]
        )
                if generated_plasmid.get_sequence() != ordered_sequence:
                    raise ValueError("Error occured, generated plasmid object does not equal sequence imported")

                self.add_plasmid(generated_plasmid)

                      
        else:
            print("Submission incomplete")

    def read_fasta(self, folder_name):
        fasta_sequences = []
        current_dict = os.getcwd()
        folder_path = os.path.join(current_dict,"storage",folder_name)

        for filename in os.listdir(folder_path):

            if filename.endswith(".fasta") or filename.endswith(".fa"):
                filepath = os.path.join(folder_path, filename)
                with open(filepath, "r") as file:
                    for record in SeqIO.parse(file, "fasta"):
                        info_pair = []
                        info_pair.append(record.id)
                        info_pair.append(str(record.seq))
                        fasta_sequences.append(info_pair)
        return fasta_sequences
    
    def read_genbank(self, folder_name):
        genbank_sequences = []
        current_dict = os.getcwd()
        folder_path = os.path.join(current_dict,"storage",folder_name)

        for filename in os.listdir(folder_path):
            if filename.endswith(".gb") or filename.endswith(".gbk"):
                filepath = os.path.join(folder_path, filename)
                with open(filepath, "r") as file:
                    for record in SeqIO.parse(file, "genbank"):
                        info_pair = []
                        info_pair.append(record.id)
                        info_pair.append(str(record.seq))
                        genbank_sequences.append(info_pair)
        return genbank_sequences
    
