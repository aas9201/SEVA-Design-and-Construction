from plasmid import Plasmid, Module, Part
from library import Library
import tkinter as tk
from tkinter import ttk
import json
from utils import import_parameters_ligation, import_insert_parameters
import textwrap
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio import SeqIO
from pydna.dseqrecord import Dseqrecord
from pydna.dseq import Dseq
from Bio import Restriction
from Bio.Restriction import RestrictionBatch
from Bio.Restriction import PacI, StyI, AvrII, BglI, SfiI, NotI, SacII, EcoRI, ApoI, Eco53kI, SacI, BanII, BanI, Acc65I, KpnI, XmaI, BsoBI, TspMI, AvaI, SmaI, BamHI, XbaI, SalI, AccI, HincII, SfcI, MlyI, PleI, SfcI, PstI, SbfI, BfuAI, BspMI, NspI, SphI, HindIII, HgaI, NotI, SpeI



def mcs_digest1(sequence1,sequence2,*enzymes):

    d_sequence_donor_plasmid = Dseq(sequence1, circular = True)
    d_sequence_recipient_plasmid = Dseq(sequence2, circular = True)

    donor_plasmid_fragments = d_sequence_donor_plasmid.cut(*enzymes)
    recipient_plasmid_fragments = d_sequence_recipient_plasmid.cut(*enzymes)

    
    d_new_sequence = Dseqrecord(donor_plasmid_fragments[0] + recipient_plasmid_fragments[1]).looped()
    new_sequence = str(d_new_sequence.seq)

    reference_point = new_sequence.find(
        "TTAATTAA"
    )  
    ordered_sequence = (
        new_sequence[reference_point:] + new_sequence[:reference_point]
    )

    digestion_data = (donor_plasmid_fragments,recipient_plasmid_fragments, ordered_sequence)

    return digestion_data



def add_flank_re(part_list,enzyme1,enzyme2):
    enzyme1_part = None
    enzyme2_part = None

    # Define the parts for each enzyme
    if enzyme1 == BamHI:
        enzyme1_part = Part("bamhi", "bamhi_id", "GGATCC", "re", "restriction enzyme site")
    elif enzyme1 == EcoRI:
        enzyme1_part = Part("ecori", "ecori_id", "GAATTC", "re", "restriction enzyme site")
    elif enzyme1 == SacI:
        enzyme1_part = Part("saci", "saci_id", "GAGCTC", "re", "restriction enzyme site")

    if enzyme2 == SalI:
        enzyme2_part = Part("sali", "sali_id", "GTCGAC", "re", "restriction enzyme site")
    elif enzyme2 == HindIII:
        enzyme2_part = Part("hindiii", "hindiii_id", "AAGCTT", "re", "restriction enzyme site")
    elif enzyme2 == PstI:
        enzyme2_part = Part("psti", "psti_id", "CTGCAG", "re", "restriction enzyme site")

    # Append the enzyme parts to the list
    if enzyme1_part:
        part_list.insert(0, enzyme1_part)  # Add to the beginning
    if enzyme2_part:
        part_list.append(enzyme2_part)  # Add to the end

    return part_list




def find_sites(plasmid_seq,enzyme1,enzyme2):
    left_site = right_site = ""

    if enzyme1 == BamHI:
        left_site = plasmid_seq.find("GGATCC")
    
    if enzyme1 == EcoRI:
        left_site = plasmid_seq.find("GAATTC")

    if enzyme1 == SacI:
        left_site = plasmid_seq.find("GAGCTC")
    
    if enzyme2 == SalI:
        right_site = [plasmid_seq.find("GTCGAC"),6]
        end = plasmid_seq.find("ACTAGT")

    
    if enzyme2 == HindIII:
        right_site = [plasmid_seq.find("AAGCTT"),6]
        end = plasmid_seq.find("ACTAGT")

    
    if enzyme2 == PstI:
        right_site = [plasmid_seq.find("CTGCAG"),6]
        end = plasmid_seq.find("ACTAGT")
    
    return [left_site,right_site, end]




class InsertConstruct:
    def __init__(self, library):
        self.library = library  
        self.root = tk.Tk()
        self.root.style = ttk.Style()
        self.root.style.theme_use('clam')
        self.root.style.theme_use('clam')  # Modern theme
        self.root.style.configure('TFrame', background='#f0f0f0')
        self.root.style.configure('TButton', font=('Helvetica', 12))
        self.root.style.configure('TLabel', background='#f0f0f0', font=('Helvetica', 12))
        self.root.style.configure('Treeview', background='white', foreground='black', fieldbackground='white')
        self.root.style.map('Treeview', background=[('selected', '#0078d7')])
        
        self.temp_selection = []
        self.plasmid_insert = []
        self.csv_data_list = []
        self.selection_value = None
        self.frame1 = ttk.Frame(self.root)
        self.frame2 = ttk.Frame(self.root)
        self.plasmid_selection_message_var = tk.StringVar()
        self.module_selection_message_var = tk.StringVar()
        self.delete_message_var = tk.StringVar()
        self.plasmid_message_var = tk.StringVar()

        for frame in (self.frame1, self.frame2):
            frame.grid(row=0, column=0, sticky='nsew')

        self.create_widgets_frame1()
        self.create_widgets_frame2()
        self.load_plasmids()
        self.load_inserts()

        self.show_frame(self.frame1)

        self.root.title("Insertion Tool")
        

        for frame in (self.frame1, self.frame2):
            frame.grid(row=0, column=0, sticky='nsew')
            self.root.grid_rowconfigure(0, weight=1)
            self.root.grid_columnconfigure(0, weight=1)

    def create_widgets_frame1(self):
        self.instructions = ttk.Label(self.frame1, text="                      Select plasmids from the list below by clicking them, then 'Select Plasmid'.\n Once selected, a new page will appear in which you can select the module you wish to insert inside.")
        self.instructions.pack(pady=10)

        self.tree = ttk.Treeview(self.frame1)
        self.tree["columns"] = ("plasmid_name", "plasmid_id", "cargo_name", "cargo_id", "marker_name", "marker_id", "replication_name", "replication_id")
        self.tree.column("#0", width=0, stretch=tk.NO)
        self.tree.column("plasmid_name", anchor=tk.W, width=120)
        self.tree.column("plasmid_id", anchor=tk.W, width=80)
        self.tree.column("cargo_name", anchor=tk.W, width=120)
        self.tree.column("cargo_id", anchor=tk.W, width=80)
        self.tree.column("marker_name", anchor=tk.W, width=120)
        self.tree.column("marker_id", anchor=tk.W, width=80)
        self.tree.column("replication_name", anchor=tk.W, width=120)
        self.tree.column("replication_id", anchor=tk.W, width=80)

        self.tree.heading("#0", text="", anchor=tk.W)
        self.tree.heading("plasmid_name", text="Plasmid Name", anchor=tk.W)
        self.tree.heading("plasmid_id", text="Plasmid ID", anchor=tk.W)
        self.tree.heading("cargo_name", text="Cargo Name", anchor=tk.W)
        self.tree.heading("cargo_id", text="Cargo ID", anchor=tk.W)
        self.tree.heading("marker_name", text="Marker Name", anchor=tk.W)
        self.tree.heading("marker_id", text="Marker ID", anchor=tk.W)
        self.tree.heading("replication_name", text="Replication Name", anchor=tk.W)
        self.tree.heading("replication_id", text="Replication ID")

        self.tree.pack(fill=tk.BOTH, expand=True)

        self.plasmid_selection_message_label = ttk.Label(self.frame1, textvariable=self.plasmid_selection_message_var)
        self.plasmid_selection_message_label.pack(pady=5)
        self.plasmid_selection_message_var.set("Please select your first plasmid, then click the 'Select' button")

        self.select_plasmid_button = ttk.Button(self.frame1, text="Select Plasmid", command=lambda: self.select_plasmid())
        self.select_plasmid_button.pack(pady=10)

        self.delete_message_label = ttk.Label(self.frame1, textvariable=self.delete_message_var)
        self.delete_message_label.pack(pady=5)

        self.delete_button = ttk.Button(self.frame1, text="Remove Last Selection", command=self.remove_plasmid_last_selected)
        self.delete_button.pack(pady=5)

        bottom_right_corner_frame = ttk.Frame(self.frame1)
        bottom_right_corner_frame.pack(side=tk.BOTTOM, anchor=tk.SE, padx=5, pady=5)
        self.display_button = ttk.Button(bottom_right_corner_frame, text="Print Added", command=self.print_selected)
        self.display_button.pack(pady=5)

        
        self.csv_frame = ttk.Frame(self.frame1)
        self.csv_frame.pack(anchor="w", padx=5, pady=5)

        self.csv_label = ttk.Label(self.csv_frame, text="CSV File Name:")
        self.csv_label.pack(side="top", anchor="nw", pady=5)

        # Pack the entry in the middle
        self.csv_entry = ttk.Entry(self.csv_frame)
        self.csv_entry.pack(side = "left", anchor="center", pady=5)

        # Pack the button at the bottom
        self.insertion_button = ttk.Button(self.csv_frame, text="Insert Plasmids", command=self.insert_plasmids)
        self.insertion_button.pack(side="bottom", anchor= "sw", pady=5)
        

    
        self.insertion_button.pack_forget()
        self.csv_entry.pack_forget()
        self.csv_label.pack_forget()



        

    def create_widgets_frame2(self):
        self.instructions2 = ttk.Label(self.frame2, text="Please select the module you wish to insert inside this plasmid, then click the 'Select' button.")
        self.instructions2.pack(pady=10)

        # TreeView for selecting inserts
        self.tree2 = ttk.Treeview(self.frame2)
        self.tree2["columns"] = ("insert_name", "insert_id", "insert_description", "part_list")
        self.tree2.column("#0", width=0, stretch=tk.NO)
        self.tree2.column("insert_name", anchor=tk.W, width=100)
        self.tree2.column("insert_id", anchor=tk.W, width=120)
        self.tree2.column("insert_description", anchor=tk.W, width=150)
        self.tree2.column("part_list", anchor=tk.W, width=500)  # Adjusted dynamically
        
        self.tree2.heading("#0", text="", anchor=tk.W)
        self.tree2.heading("insert_name", text="Inert Name", anchor=tk.W)
        self.tree2.heading("insert_id", text="Insert ID", anchor=tk.W)
        self.tree2.heading("insert_description", text="Insert Description", anchor=tk.W)
        self.tree2.heading("part_list", text="Part List", anchor=tk.W)


        self.tree2.pack(fill=tk.BOTH, expand=True)

        # Buttons at the bottom
        options = ["BamHI","EcoRI","SacI"]
        self.selection = ttk.Combobox(self.frame2, values=options, state="readonly", width=45)
        self.selection.set("Select Left Flanking Enzyme Recognition Site")
        self.selection.pack(pady=5)

        options2 = ["PstI","SalI","HindIII"]
        self.selection2 = ttk.Combobox(self.frame2, values=options2, state="readonly", width=45)
        self.selection2.set("Select Right Flanking Enzyme Recognition Site")
        self.selection2.pack(pady=5)

        # Input Frame for Plasmid Information
        self.input_frame = ttk.Frame(self.frame2)
        self.input_frame.pack(pady=5)  # Ensure it appears below the options

        ttk.Label(self.input_frame, text="Plasmid Name:").grid(row=0, column=0, padx=5, pady=5)
        self.plasmid_name_entry = ttk.Entry(self.input_frame)
        self.plasmid_name_entry.grid(row=0, column=1, padx=5, pady=5)

        ttk.Label(self.input_frame, text="Plasmid ID:").grid(row=1, column=0, padx=5, pady=5)
        self.plasmid_id_entry = ttk.Entry(self.input_frame)
        self.plasmid_id_entry.grid(row=1, column=1, padx=5, pady=5)

        ttk.Label(self.input_frame, text="Plasmid Description:").grid(row=2, column=0, padx=5, pady=5)
        self.plasmid_description_entry = ttk.Entry(self.input_frame)
        self.plasmid_description_entry.grid(row=2, column=1, padx=5, pady=5)


        self.select_insert_button = ttk.Button(self.frame2, text="Select Insert", command=lambda: self.select_insert())
        self.select_insert_button.pack(pady=10)

        self.back_button = ttk.Button(self.frame2, text="Go back and reselect plasmid", command=lambda: self.go_back())
        self.back_button.pack(side="bottom", pady=10)
        

    def show_frame(self, frame):
        frame.tkraise()

    def select_plasmid(self):

        selected_items = self.tree.selection() 
        unique_id  = str(selected_items[0])
        plasmid =  self.library.retrieve_plasmid(unique_id)
        self.temp_selection.append(plasmid)
        self.show_frame(self.frame2)
        self.module_selection_message_var.set(textwrap.dedent(f"Please select the module you wish to insert inside the plasmid: '{plasmid.name}', '{plasmid.unique_id}'"))
        self.delete_message_var.set("")

    def select_insert(self):
        enzyme_sequences = {
                "BamHI": "GGATCC",
                "EcoRI": "GAATTC",
                "SacI": "GAGCTC",
                "PstI": "CTGCAG",
                "SalI": "GTCGAC",
                "HindIII": "AAGCTT"
            }
        selected_insert = self.tree2.selection()
        if selected_insert:  # Ensure there is at least one selected item
            insert_id = selected_insert[0]  
            insert_data = self.library.retrieve_insert(insert_id)
            print(insert_data)
            enzyme1 = self.selection.get()
            enzyme2 = self.selection2.get()

        part_list = []
        for i in insert_data['structure']:
            part = self.library.retrieve_part(i)
            part_list.append(part)
        
        insert_seq = ""
        for i in part_list:
            insert_seq += i.sequence
        
        site1 = enzyme_sequences.get(enzyme1, "")
        site2 = enzyme_sequences.get(enzyme2, "")

        site1_present = site1 in insert_seq
        site2_present = site2 in insert_seq

        if site1_present:
            print(f"Insert contains forbidden sites: {enzyme1}. Please choose a different insert or modify.")
            return 

        if site2_present:
            print(f"Insert contains forbidden sites: {enzyme2}. Please choose a different insert or modify.")
            return
        
        existing_plasmid = self.library._plasmid_exists(self.plasmid_id_entry.get())
        if existing_plasmid:
            print(f"Plasmid with Unique ID: {self.plasmid_id_entry.get()} is already within database.")
        else:
            self.temp_selection.append(insert_data)
            self.temp_selection.append(self.plasmid_name_entry.get())
            self.temp_selection.append(self.plasmid_id_entry.get())
            self.temp_selection.append(self.plasmid_description_entry.get())
            self.temp_selection.append(enzyme1)
            self.temp_selection.append(enzyme2)
            self.plasmid_insert.append(self.temp_selection)
            self.temp_selection = []
            self.show_frame(self.frame1)
            print(self.plasmid_insert)
            part_list = []
            enzyme1 = None
            enzyme2 = None
            self.show_frame(self.frame1)
            self.clear_input_fields()
            self.csv_label.pack(pady=5)
            self.csv_entry.pack(pady = 5)
            self.insertion_button.pack(pady=5)
            self.delete_message_var.set("")
            self.temp_selection = []




    def clear_input_fields(self):
        self.plasmid_name_entry.delete(0, tk.END)
        self.plasmid_id_entry.delete(0, tk.END)
        self.plasmid_description_entry.delete(0, tk.END)

    
    def select(self, event):
        self.selection_value = self.selection.get()
        print(self.selection_value)



    def insert_plasmids(self):
        csv_name = self.csv_entry.get()
        for pairs in self.plasmid_insert:
            construct_info = []
            constructed_plasmid = self.create_new_module(pairs[1], pairs[0], pairs[2], pairs[3], pairs[4], eval(pairs[5]), eval(pairs[6]))
            print(constructed_plasmid[0])
            construct_info.append(pairs[1]['name'])
            reaction_frag_donor = constructed_plasmid[2][0]
            reaction_frag_recipient = constructed_plasmid[1][0]
            discarded_frag_recipient = constructed_plasmid[1][1]
            construct_info = [pairs[1]['name'],len(reaction_frag_donor),constructed_plasmid[0].name,len(reaction_frag_recipient),len(discarded_frag_recipient),pairs[5],pairs[6]]
            self.csv_data_list.append(construct_info)

        
        import_insert_parameters(csv_name,self.csv_data_list)
        import_parameters_ligation(csv_name,self.csv_data_list,"module")





    def load_plasmids(self):
        data = self.library._retrieve_all_json_plasmids()
        for unique_id, plasmid_json in data:
            plasmid_data = json.loads(plasmid_json)
            plasmid_id = plasmid_data["unique_id"]
            plasmid_name = plasmid_data["name"]
            
            cargo_id = cargo_name = None
            marker_id = marker_name = None
            replication_id = replication_name = None

            has_desired_cargo = False
            
            for module in plasmid_data["structure"]:
                if module["module_type"] == "cargo":
                    cargo_id = module["unique_id"]
                    cargo_name = module["name"]
                    if cargo_id == "mcs":
                        has_desired_cargo = True
                elif module["module_type"] == "marker":
                    marker_id = module["unique_id"]
                    marker_name = module["name"]
                elif module["module_type"] == "replication":
                    replication_id = module["unique_id"]
                    replication_name = module["name"]

            if has_desired_cargo:
                if self.tree.exists(plasmid_id):
                    self.tree.item(plasmid_id, values=(plasmid_name, plasmid_id, cargo_name, cargo_id, marker_name, marker_id, replication_name, replication_id))
                else:
                    self.tree.insert("", "end", iid=plasmid_id, values=(plasmid_name, plasmid_id, cargo_name, cargo_id, marker_name, marker_id, replication_name, replication_id))


    def load_inserts(self):
        data = self.library._retrieve_all_json_inserts()
        part_list = insert_id = insert_name = insert_description = ""

        for unique_id, insert_json in data:
            insert_data = json.loads(insert_json)

            insert_id = insert_data["unique_id"]
            insert_name = insert_data["name"]
            insert_description = insert_data["description"]
            part_list = insert_data["structure"]


        
            max_part_list_width = max(250, len(part_list) * 7) 
            self.tree2.column("part_list", width=max_part_list_width)
            
            if self.tree2.exists(insert_id):
                self.tree2.item(insert_id, values=(insert_name, insert_id, insert_description, part_list))
            else:
                self.tree2.insert("", "end", iid=insert_id, values=(insert_name, insert_id, insert_description, part_list))


    def remove_plasmid_last_selected(self):
        if self.plasmid_insert:
            self.delete_message_var.set(f"Plasmid: {self.plasmid_insert[-1][0].name} and {self.plasmid_insert[-1][1]['name']} have been removed")
            self.plasmid_insert.pop()
            

    def go_back(self):
        if self.temp_selection:
            self.temp_selection.pop()
            self.show_frame(self.frame1)

    def print_selected(self):
        print("The Plasmid and Insert pairs which have been selected:")
        for i in self.plasmid_insert:
            print(f"Plasmid: {i[0].name}, Insert: {i[1]['name']}")

    def create_new_module(self,insert_dict,plasmid,new_plasmid_name, new_plasmid_id, new_plasmid_description, *enzymes):
        part_list = []
        for i in insert_dict['structure']:
            part = self.library.retrieve_part(i)
            part_list.append(part)

        insert_list = add_flank_re(part_list,*enzymes)
        insert_seq = ""
        for i in insert_list:
            insert_seq += i.sequence

        new_plasmid_model = mcs_digest1(insert_seq,plasmid.get_sequence(),*enzymes)
        new_plasmid_seq = new_plasmid_model[2]
        new_plasmid_model_donor_fragments = new_plasmid_model[0]
        new_plasmid_model_recipient_fragments = new_plasmid_model[1]
        
        print("below is  seq")
        print(new_plasmid_seq)
        sites = find_sites(new_plasmid_seq,*enzymes)

        misc_part1seq = new_plasmid_seq[8:sites[0]]
        misc_part2seq = new_plasmid_seq[sites[1][0]+sites[1][1]:sites[2]]

        misc_part1 = Part(f"mcs_left_{plasmid.name}", f"mcs_left_{plasmid.unique_id}", misc_part1seq, "misc", "sequence left after insertion")
        misc_part2 = Part(f"mcs_right_{plasmid.name}", f"mcs_right_{plasmid.unique_id}", misc_part2seq, "misc", "sequence left after insertion")
    
        insert_list.insert(0,misc_part1)
        insert_list.append(misc_part2)


        new_module = Module(f"{insert_dict['name']}", f"{insert_dict['unique_id']}","cargo", f"{insert_dict['description']}", structure=insert_list)

        new_plasmid = Plasmid(f"{new_plasmid_name}",f"{new_plasmid_id}", f"{new_plasmid_description}")
        new_plasmid.structure = [new_module,plasmid.structure[2],plasmid.structure[4]]

        if plasmid.scar_info:
            scar_data = plasmid.scar_info
            for scar_id, scar_info in scar_data.items():
                scar_part = Part(
                    name=scar_info['name'],
                    unique_id=scar_id,
                    sequence=scar_info['sequence'],
                    role='scar',
                    description=scar_info['description']
                )
                new_plasmid.insert_scar(scar_part, scar_info['module'], scar_info['position'])

        return new_plasmid, new_plasmid_model_recipient_fragments, new_plasmid_model_donor_fragments

    def run(self):
        self.root.mainloop()



def insertion_tool(library):  
    builder = InsertConstruct(library) 
    builder.run()  



from plasmid import Plasmid, Module, Part
from library import Library
import tkinter as tk
from tkinter import ttk
import json
from utils import process_plasmid_list, import_single_parameters, import_module_single_parameters, import_module_multi_parameters, import_multi_parameters, process_module_plasmid_list, import_parameters_ligation
import textwrap


class ModuleConstruct:
    def __init__(self, library, module_type):
        self.library = library  
        self.root = tk.Tk()
        self.root.style = ttk.Style()
        self.root.style.theme_use('clam')
        self.root.style.theme_use('clam')  # Modern theme
        self.root.style.configure('TFrame', background='#f0f0f0')
        self.root.style.configure('TButton', font=('Helvetica', 12))
        self.root.style.configure('TLabel', background='#f0f0f0', font=('Helvetica', 12))
        self.root.style.configure('Treeview', background='white', foreground='black', fieldbackground='white')
        self.root.style.map('Treeview', background=[('selected', '#0078d7')])
        
        self.temp_selection = []
        self.plasmid_module = []
        self.selection_value = None
        self.frame1 = ttk.Frame(self.root)
        self.frame2 = ttk.Frame(self.root)
        self.plasmid_selection_message_var = tk.StringVar()
        self.module_selection_message_var = tk.StringVar()
        self.delete_message_var = tk.StringVar()
        self.plasmid_message_var = tk.StringVar()
        self.module_type = module_type
        for frame in (self.frame1, self.frame2):
            frame.grid(row=0, column=0, sticky='nsew')

        self.create_widgets_frame1()
        self.create_widgets_frame2()
        self.load_plasmids()
        self.load_modules()

        self.show_frame(self.frame1)

        if module_type == "cargo":
            self.root.title("Cargo Module Insertion")
        elif module_type == "marker":
            self.root.title("Marker Module Insertion")
        elif module_type == "replication":
            self.root.title("Replication Module Insertion")
        else:
            raise ValueError("Please ensure module type is either 'cargo', 'marker' or 'replication'")
        

        for frame in (self.frame1, self.frame2):
            frame.grid(row=0, column=0, sticky='nsew')
            self.root.grid_rowconfigure(0, weight=1)
            self.root.grid_columnconfigure(0, weight=1)

    def create_widgets_frame1(self):
        self.instructions = ttk.Label(self.frame1, text="                      Select plasmids from the list below by clicking them, then 'Select Plasmid'.\n Once selected, a new page will appear in which you can select the module you wish to insert inside.")
        self.instructions.pack(pady=10)

        self.tree = ttk.Treeview(self.frame1)
        self.tree["columns"] = ("plasmid_name", "plasmid_id", "cargo_name", "cargo_id", "marker_name", "marker_id", "replication_name", "replication_id")
        self.tree.column("#0", width=0, stretch=tk.NO)
        self.tree.column("plasmid_name", anchor=tk.W, width=120)
        self.tree.column("plasmid_id", anchor=tk.W, width=80)
        self.tree.column("cargo_name", anchor=tk.W, width=120)
        self.tree.column("cargo_id", anchor=tk.W, width=80)
        self.tree.column("marker_name", anchor=tk.W, width=120)
        self.tree.column("marker_id", anchor=tk.W, width=80)
        self.tree.column("replication_name", anchor=tk.W, width=120)
        self.tree.column("replication_id", anchor=tk.W, width=80)

        self.tree.heading("#0", text="", anchor=tk.W)
        self.tree.heading("plasmid_name", text="Plasmid Name", anchor=tk.W)
        self.tree.heading("plasmid_id", text="Plasmid ID", anchor=tk.W)
        self.tree.heading("cargo_name", text="Cargo Name", anchor=tk.W)
        self.tree.heading("cargo_id", text="Cargo ID", anchor=tk.W)
        self.tree.heading("marker_name", text="Marker Name", anchor=tk.W)
        self.tree.heading("marker_id", text="Marker ID", anchor=tk.W)
        self.tree.heading("replication_name", text="Replication Name", anchor=tk.W)
        self.tree.heading("replication_id", text="Replication ID")

        self.tree.pack(fill=tk.BOTH, expand=True)

        self.plasmid_selection_message_label = ttk.Label(self.frame1, textvariable=self.plasmid_selection_message_var)
        self.plasmid_selection_message_label.pack(pady=5)
        self.plasmid_selection_message_var.set("Please select your first plasmid, then click the 'Select' button")

        self.select_plasmid_button = ttk.Button(self.frame1, text="Select Plasmid", command=lambda: self.select_plasmid())
        self.select_plasmid_button.pack(pady=10)

        self.delete_message_label = ttk.Label(self.frame1, textvariable=self.delete_message_var)
        self.delete_message_label.pack(pady=5)

        self.delete_button = ttk.Button(self.frame1, text="Remove Last Selection", command=self.remove_plasmid_last_selected)
        self.delete_button.pack(pady=5)

        bottom_right_corner_frame = ttk.Frame(self.frame1)
        bottom_right_corner_frame.pack(side=tk.BOTTOM, anchor=tk.SE, padx=5, pady=5)
        self.display_button = ttk.Button(bottom_right_corner_frame, text="Print Added", command=self.print_selected)
        self.display_button.pack(pady=5)

        self.input_frame = ttk.Frame(self.frame1)
        self.input_frame.pack(pady=5)

        self.plasmid_message_label = ttk.Label(self.frame1, textvariable=self.plasmid_message_var)
        self.plasmid_message_label.pack(pady=5)

        ttk.Label(self.input_frame, text="Plasmid Name:").grid(row=0, column=0, padx=5, pady=5)
        self.plasmid_name_entry = ttk.Entry(self.input_frame)
        self.plasmid_name_entry.grid(row=0, column=1, padx=5, pady=5)

        ttk.Label(self.input_frame, text="Plasmid ID:").grid(row=1, column=0, padx=5, pady=5)
        self.plasmid_id_entry = ttk.Entry(self.input_frame)
        self.plasmid_id_entry.grid(row=1, column=1, padx=5, pady=5)

        ttk.Label(self.input_frame, text="Plasmid Description:").grid(row=2, column=0, padx=5, pady=5)
        self.plasmid_description_entry = ttk.Entry(self.input_frame)
        self.plasmid_description_entry.grid(row=2, column=1, padx=5, pady=5)
    
        self.submit_button = ttk.Button(self.input_frame, text="Submit", command=self.submit_plasmid)
        self.submit_button.grid(row=3, column=0, columnspan=2, pady=5)
        
        self.csv_frame = ttk.Frame(self.frame1)
        self.csv_frame.pack(anchor="w", padx=5, pady=5)

        self.csv_label = ttk.Label(self.csv_frame, text="CSV File Name:")
        self.csv_label.pack(side="top", anchor="nw", pady=5)

        # Pack the entry in the middle
        self.csv_entry = ttk.Entry(self.csv_frame)
        self.csv_entry.pack(side = "left", anchor="center", pady=5)

        # Pack the button at the bottom
        self.insertion_button = ttk.Button(self.csv_frame, text="Insert Plasmids", command=self.insert_plasmids)
        self.insertion_button.pack(side="bottom", anchor= "sw", pady=5)
        
        options = ["Single Channel", "Multi-Channel"]
        self.selection = ttk.Combobox(self.frame1, values=options, state= "readonly")
        self.selection.set("Select Protocol Type")
        self.selection.bind("<<ComboboxSelected>>", self.select)
        self.selection.pack(pady=5)

        self.input_frame.pack_forget()
        self.insertion_button.pack_forget()
        self.csv_entry.pack_forget()
        self.csv_label.pack_forget()



        

    def create_widgets_frame2(self):
        self.instructions2 = ttk.Label(self.frame2, text="Please select the module you wish to insert inside this plasmid, then click the 'Select' button.")
        self.instructions2.pack(pady=10)

        self.tree2 = ttk.Treeview(self.frame2)
        self.tree2["columns"] = ("module_name", "module_id", "module_description", "part_list")
        self.tree2.column("#0", width=0, stretch=tk.NO)
        self.tree2.column("module_name", anchor=tk.W, width=100)
        self.tree2.column("module_id", anchor=tk.W, width=120)
        self.tree2.column("module_description", anchor=tk.W, width=150)
        self.tree2.column("part_list", anchor=tk.W, width=500)  # Initial width, adjust based on content

        self.tree2.heading("#0", text="", anchor=tk.W)

        if self.module_type == "cargo":
            self.tree2.heading("module_name", text="Cargo Name", anchor=tk.W)
            self.tree2.heading("module_id", text="Cargo ID", anchor=tk.W)
            self.tree2.heading("module_description", text="Cargo Description", anchor=tk.W)
            self.tree2.heading("part_list", text="Part List", anchor=tk.W)

        if self.module_type == "marker":
            self.tree2.heading("module_name", text="Marker Name", anchor=tk.W)
            self.tree2.heading("module_id", text="Marker ID", anchor=tk.W)
            self.tree2.heading("module_description", text="Marker Description", anchor=tk.W)
            self.tree2.heading("part_list", text="Part List", anchor=tk.W)

        if self.module_type == "replication":
            self.tree2.heading("module_name", text="Replication Name", anchor=tk.W)
            self.tree2.heading("module_id", text="Repilcation ID", anchor=tk.W)
            self.tree2.heading("module_description", text="Replication Description", anchor=tk.W)
            self.tree2.heading("part_list", text="Part List", anchor=tk.W)


        self.tree2.pack(fill=tk.BOTH, expand=True)

        

        self.module_selection_message_label = ttk.Label(self.frame2, textvariable=self.module_selection_message_var)
        self.module_selection_message_label.pack(pady=5)

        self.back_button = ttk.Button(self.frame2, text="Go back and reselect plasmid", command=lambda: self.go_back())
        self.back_button.pack(side = "bottom", anchor = "s" ,pady=10)

        self.select_module_button = ttk.Button(self.frame2, text="Select Module", command=lambda: self.select_module())
        self.select_module_button.pack(side = "bottom", anchor = "s" ,pady=10)

        

    def show_frame(self, frame):
        frame.tkraise()

    def select_plasmid(self):

        selected_items = self.tree.selection() 
        unique_id  = str(selected_items[0])
        plasmid =  self.library.retrieve_plasmid(unique_id)
        self.temp_selection.append(plasmid)
        self.show_frame(self.frame2)
        self.module_selection_message_var.set(textwrap.dedent(f"Please select the module you wish to insert inside the plasmid: '{plasmid.name}', '{plasmid.unique_id}'"))
        self.delete_message_var.set("")

    def select_module(self):

        selected_items = self.tree2.selection() 
        unique_id  = str(selected_items[0])
        module =  self.library.retrieve_module(unique_id)
        self.temp_selection.append(module)
        self.plasmid_module.append(self.temp_selection)
        print(f"Plasmid: {self.plasmid_module[-1][0].name} and Marker: {self.plasmid_module[-1][1].name} have been added")
        if self.plasmid_module:
            self.plasmid_selection_message_var.set(           "Please select your next plasmid and module pair\n If you are finished, name your file and click the insert button")



        if self.module_type == "cargo":
            if len(self.temp_selection) == 2:
                module_id_list = [self.temp_selection[1].unique_id,self.temp_selection[0].structure[2].unique_id,self.temp_selection[0].structure[4].unique_id]
                module_id_check = self.library._retrieve_all_module_id()
                self.csv_label.pack(pady=5)
                self.csv_entry.pack(pady = 5)
                self.insertion_button.pack(pady=5)

        if self.module_type == "marker":
            if len(self.temp_selection) == 2:
                module_id_list = [self.temp_selection[0].structure[0],self.temp_selection[1].unique_id,self.temp_selection[0].structure[4].unique_id]
                module_id_check = self.library._retrieve_all_module_id()
                self.csv_label.pack(pady=5)
                self.csv_entry.pack(pady = 5)
                self.insertion_button.pack(pady=5)

        if self.module_type == "replication":
            if len(self.temp_selection) == 2:
                module_id_list = [self.temp_selection[0].structure[0],self.temp_selection[0].structure[2].unique_id,self.temp_selection[1].unique_id]
                module_id_check = self.library._retrieve_all_module_id()
                self.csv_label.pack(pady=5)
                self.csv_entry.pack(pady = 5)
                self.insertion_button.pack(pady=5)
                


        if module_id_list not in module_id_check:
                self.plasmid_selection_message_var.set("")
                self.input_frame.pack(pady=5)
                self.select_plasmid_button.pack_forget()
                self.delete_button.pack_forget()
                self.selection.pack_forget()
                self.plasmid_message_var.set(f"Generated Plasmid is new, please insert details below to add to database")
                self.insertion_button.pack_forget()
                self.csv_entry.pack_forget()
                self.csv_label.pack_forget()

    
        self.show_frame(self.frame1)
        self.delete_message_var.set("")
        self.temp_selection = []




    def submit_plasmid(self):
        plasmid_name = self.plasmid_name_entry.get()
        plasmid_id = self.plasmid_id_entry.get()
        plasmid_description = self.plasmid_description_entry.get()
        print(self.plasmid_module[-1][1], self.plasmid_module[-1][0].structure[2],self.plasmid_module[-1][0].structure[4])

        if self.module_type == "cargo":
            plasmid = Plasmid(plasmid_name, plasmid_id, plasmid_description)
            plasmid.structure = [self.plasmid_module[-1][1], self.plasmid_module[-1][0].structure[2],self.plasmid_module[-1][0].structure[4]]

        if self.module_type == "marker":
            plasmid = Plasmid(plasmid_name, plasmid_id, plasmid_description)
            plasmid.structure = [self.plasmid_module[-1][0].structure[0], self.plasmid_module[-1][1],self.plasmid_module[-1][0].structure[4]]
        
        if self.module_type == "replication":
            plasmid = Plasmid(plasmid_name, plasmid_id, plasmid_description)
            plasmid.structure = [self.plasmid_module[-1][0].structure[0], self.plasmid_module[-1][0].structure[2],self.plasmid_module[-1][1]]
    

        if self.plasmid_module[-1][0].scar_info:
            scar_data = self.plasmid_objects[-1][1].scar_info
            for scar_id, scar_info in scar_data.items():
                scar_part = Part(
                    name=scar_info['name'],
                    unique_id=scar_id,
                    sequence=scar_info['sequence'],
                    role='scar',
                    description=scar_info['description']
                )
                plasmid.insert_scar(scar_part, scar_info['module'], scar_info['position'])

        checker = self.library._constructor_add_plasmid(plasmid)

        if checker:
            self.input_frame.pack_forget()
            self.select_plasmid_button.pack(pady=5)
            self.delete_button.pack(pady=5)
            self.plasmid_message_var.set("")
            self.csv_label.pack(pady=5)
            self.csv_entry.pack(pady = 5)
            self.insertion_button.pack(pady=5)
            self.clear_input_fields()
            self.selection.pack(pady=5)
            self.plasmid_selection_message_var.set("           Please select your next plasmid and module pair\n If you are finished, name your file and click the insert button")

    def clear_input_fields(self):
        self.plasmid_name_entry.delete(0, tk.END)
        self.plasmid_id_entry.delete(0, tk.END)
        self.plasmid_description_entry.delete(0, tk.END)

    
    def select(self, event):
        self.selection_value = self.selection.get()
        print(self.selection_value)



    def insert_plasmids(self):
        csv_name = self.csv_entry.get()
        module_to_plasmid = self.library._module_to_plasmid()

        if self.module_type == "cargo":
            if self.selection_value: #Just to check if single/multi-channel chosen
                for pairs in self.plasmid_module: #Just do add plasmid_id to the plasmid_object list for reference
                    recipient_plasmid = pairs[0]
                    module_inserting = pairs[1]
                    desired_plasmid_module_list = (module_inserting.unique_id,recipient_plasmid.structure[2].unique_id,recipient_plasmid.structure[4].unique_id)
                    desired_plasmid_id = module_to_plasmid[desired_plasmid_module_list]
                    pairs.append(desired_plasmid_id)

        if self.module_type == "marker":
            if self.selection_value: #Just to check if single/multi-channel chosen
                for pairs in self.plasmid_module: #Just do add plasmid_id to the plasmid_object list for reference
                    recipient_plasmid = pairs[0]
                    module_inserting = pairs[1]
                    desired_plasmid_module_list = (recipient_plasmid.structure[0].unique_id,module_inserting.unique_id,recipient_plasmid.structure[4].unique_id)
                    desired_plasmid_id = module_to_plasmid[desired_plasmid_module_list]
                    pairs.append(desired_plasmid_id)
            
        if self.module_type == "replication":
            if self.selection_value: #Just to check if single/multi-channel chosen
                for pairs in self.plasmid_module: #Just do add plasmid_id to the plasmid_object list for reference
                    recipient_plasmid = pairs[0]
                    module_inserting = pairs[1]
                    desired_plasmid_module_list = (recipient_plasmid.structure[0].unique_id,recipient_plasmid.structure[2].unique_id,module_inserting.unique_id)
                    desired_plasmid_id = module_to_plasmid[desired_plasmid_module_list]
                    pairs.append(desired_plasmid_id)






        if self.selection_value == "Single Channel":
            input_csv = process_module_plasmid_list(self.plasmid_module,self.module_type)  
            import_module_single_parameters(csv_name,input_csv)
            import_parameters_ligation(csv_name,input_csv,"module")
        elif self.selection_value =="Multi-Channel":
            input_csv = process_module_plasmid_list(self.plasmid_module,self.module_type)
            import_module_multi_parameters(csv_name,input_csv)
            import_parameters_ligation(csv_name,input_csv,"module")
        elif self.selection_value == None:
            print("Please select a protocol type")




    def load_plasmids(self):
        data = self.library._retrieve_all_json_plasmids()
        for unique_id, plasmid_json in data:
            plasmid_data = json.loads(plasmid_json)
            plasmid_id = plasmid_data["unique_id"]
            plasmid_name = plasmid_data["name"]
            
            for module in plasmid_data["structure"]:
                if module["module_type"] == "cargo":
                    cargo_id = module["unique_id"]
                    cargo_name = module["name"]
                if module["module_type"] == "marker":
                    marker_id = module["unique_id"]
                    marker_name = module["name"]
                if module["module_type"] == "replication":
                    replication_id = module["unique_id"]
                    replication_name = module["name"]

            if self.tree.exists(plasmid_id):
                self.tree.item(plasmid_id, values=(plasmid_name, plasmid_id, cargo_name, cargo_id, marker_name, marker_id, replication_name, replication_id))
            else:
                self.tree.insert("", "end", iid=plasmid_id, values=(plasmid_name, plasmid_id, cargo_name, cargo_id, marker_name, marker_id, replication_name, replication_id))

    def load_modules(self):
        data = self.library._retrieve_all_json_modules()
        part_list_str = module_id = module_name = module_description = ""

        for unique_id, module_json in data:
            module_data = json.loads(module_json)

            if module_data["module_type"] == self.module_type:
                module_id = module_data["unique_id"]
                module_name = module_data["name"]
                module_description = module_data["description"]
                part_list = [part["name"] for part in module_data["structure"]]
                part_list_str = ", ".join(part_list)

            

            # Dynamically adjust the column width based on content
            max_part_list_width = max(250, len(part_list_str) * 7)  # 8 pixels per character
            self.tree2.column("part_list", width=max_part_list_width)
            
            if self.tree2.exists(module_id):
                self.tree2.item(module_id, values=(module_name, module_id, module_description, part_list_str))
            else:
                self.tree2.insert("", "end", iid=module_id, values=(module_name, module_id, module_description, part_list_str))


    def remove_plasmid_last_selected(self):
        if self.plasmid_module:
            self.delete_message_var.set(f"Plasmid: {self.plasmid_module[-1][0].name} and {self.plasmid_module[-1][1].name} have been removed")
            self.plasmid_module.pop()
            

    def go_back(self):
        if self.temp_selection:
            self.temp_selection.pop()
            self.show_frame(self.frame1)

    def print_selected(self):
        print("The Plasmids which have been selected:")
        for i in self.plasmid_module:
            print(f"Plasmid: {i[0].name}, Module: {i[1].name}")



    def run(self):
        self.root.mainloop()



def module_insertion_tool(library, module_type):  
    builder = ModuleConstruct(library, module_type) 
    builder.run()  



from plasmid import Plasmid, Module, Part
import tkinter as tk
from tkinter import ttk
import json
from utils import process_plasmid_list, import_single_parameters, import_multi_parameters, import_parameters_ligation

class PlasmidConstruct:
    def __init__(self, library, module_type):
        self.library = library 
        self.module_type = module_type
        self.root = tk.Tk() 
        self.root.style = ttk.Style()
        self.root.style.theme_use('clam')
        self.root.style.configure('TFrame', background='#f0f0f0')
        self.root.style.configure('TButton', font=('Helvetica', 12))
        self.root.style.configure('TLabel', background='#f0f0f0', font=('Helvetica', 12))
        self.root.style.configure('Treeview', background='white', foreground='black', fieldbackground='white')
        self.root.style.map('Treeview', background=[('selected', '#0078d7')])
        self.message_var = tk.StringVar() 
        self.message_var.set("Please select the first Donor Plasmid ")
        self.delete_message_var = tk.StringVar() 
        self.plasmid_message_var = tk.StringVar()
        self.create_widgets()
        self.load_plasmids()
        self.plasmid_objects = []  
        self.temp_selection = []
        self.selection_value = None
    

        if module_type == "cargo":
            self.root.title("Cargo Module Insertion")
        elif module_type == "marker":
            self.root.title("Marker Module Insertion")
        elif module_type == "replication":
            self.root.title("Replication Module Insertion")
        else:
            raise ValueError("Please ensure module type is either 'cargo', 'marker' or 'replication'")

        self.root.grid_rowconfigure(0, weight=1)
        self.root.grid_columnconfigure(0, weight=1)

    def create_widgets(self):
        self.main_frame = ttk.Frame(self.root, padding="10")
        self.main_frame.grid(row=0, column=0, sticky="nsew")

        if self.module_type =="cargo":
            self.instructions = ttk.Label(self.main_frame, text="  Select plasmids from the list below by clicking them, then 'Select Plasmid'.\n      Click Donor and Recipient plasmid pairs until you are then finished,\n           then click the 'Insert Plasmids' button to generate the CSV file.\n\n                               Donor: Where Cargo module is going from.\n                               Recipient: Where Cargo module is going to")
            self.instructions.pack(pady=10)
        if self.module_type =="marker":
            self.instructions = ttk.Label(self.main_frame, text="  Select plasmids from the list below by clicking them, then 'Select Plasmid'.\n      Click Donor and Recipient plasmid pairs until you are then finished,\n           then click the 'Insert Plasmids' button to generate the CSV file.\n\n                               Donor: Where Marker module is going from.\n                               Recipient: Where Marker module is going to")
            self.instructions.pack(pady=10)
        if self.module_type =="replication":
            self.instructions = ttk.Label(self.main_frame, text="  Select plasmids from the list below by clicking them, then 'Select Plasmid'.\n      Click Donor and Recipient plasmid pairs until you are then finished,\n           then click the 'Insert Plasmids' button to generate the CSV file.\n\n                               Donor: Where Replication module is going from.\n                               Recipient: Where Replication module is going to")
            self.instructions.pack(pady=10)

        tree_frame = ttk.Frame(self.main_frame)
        tree_frame.pack(fill=tk.BOTH, expand=True)

        tree_scroll = ttk.Scrollbar(tree_frame)
        tree_scroll.pack(side=tk.RIGHT, fill=tk.Y)

        self.tree = ttk.Treeview(tree_frame, yscrollcommand=tree_scroll.set)
        tree_scroll.config(command=self.tree.yview)
        self.tree["columns"] = ("plasmid_name", "plasmid_id", "cargo_name", "cargo_id", "marker_name", "marker_id", "replication_name", "replication_id")
        self.tree.column("#0", width=0, stretch=tk.NO)
        self.tree.column("plasmid_name", anchor=tk.W, width=120)
        self.tree.column("plasmid_id", anchor=tk.W, width=80)
        self.tree.column("cargo_name", anchor=tk.W, width=120)
        self.tree.column("cargo_id", anchor=tk.W, width=80)
        self.tree.column("marker_name", anchor=tk.W, width=120)
        self.tree.column("marker_id", anchor=tk.W, width=80)
        self.tree.column("replication_name", anchor=tk.W, width=120)
        self.tree.column("replication_id", anchor=tk.W, width=80)

        self.tree.heading("#0", text="", anchor=tk.W)
        self.tree.heading("plasmid_name", text=" Plasmid Name", anchor=tk.W)
        self.tree.heading("plasmid_id", text="Plasmid ID", anchor=tk.W)
        self.tree.heading("cargo_name", text=" Cargo Name", anchor=tk.W)
        self.tree.heading("cargo_id", text="Cargo ID", anchor=tk.W)
        self.tree.heading("marker_name", text=" Marker Name", anchor=tk.W)
        self.tree.heading("marker_id", text="Marker ID", anchor=tk.W)
        self.tree.heading("replication_name", text=" Replication Name", anchor=tk.W)
        self.tree.heading("replication_id", text="Replication ID", anchor=tk.W)

        self.tree.pack(fill=tk.BOTH, expand=True)

        self.message_label = ttk.Label(self.main_frame, textvariable=self.message_var)
        self.message_label.pack(pady=5)

        self.delete_message_label = ttk.Label(self.main_frame, textvariable=self.delete_message_var)
        self.delete_message_label.pack(pady=5)

        self.select_button = ttk.Button(self.main_frame, text="Select Plasmid", command=self.select_plasmid)
        self.select_button.pack(pady=5)

        self.delete_button = ttk.Button(self.main_frame, text="Delete Last Selection", command=self.remove_last_selected)
        self.delete_button.pack(pady=5)

        bottom_right_corner_frame = ttk.Frame(self.main_frame)
        bottom_right_corner_frame.pack(side=tk.BOTTOM, anchor=tk.SE, padx=5, pady=5)
        self.display_button = ttk.Button(bottom_right_corner_frame, text="Print Added Plasmids", command=self.print_selected)
        self.display_button.pack(pady=5)

        self.selection_frame = ttk.Frame(self.main_frame)
        self.selection_frame.pack(anchor="ne", padx=5, pady=5)
        options = ["Single Channel", "Multi-Channel"]
        self.selection = ttk.Combobox(self.selection_frame, values=options, state= "readonly")
        self.selection.set("Select Protocol Type")
        self.selection.bind("<<ComboboxSelected>>", self.select)
        self.selection.pack(pady=5)

        self.csv_frame = ttk.Frame(self.main_frame)
        self.csv_frame.pack(anchor="w", padx=5, pady=5)

        # Pack the label at the top
        self.csv_label = ttk.Label(self.csv_frame, text="CSV File Name:")
        self.csv_label.pack(side="top", anchor="nw", pady=5)

        # Pack the entry in the middle
        self.csv_entry = ttk.Entry(self.csv_frame)
        self.csv_entry.pack(side = "left", anchor="center", pady=5)

        # Pack the button at the bottom
        self.insertion_button = ttk.Button(self.csv_frame, text="Insert Plasmids", command=self.insert_plasmids)
        self.insertion_button.pack(side="bottom", anchor= "sw", pady=5)

        self.plasmid_message_label = ttk.Label(self.main_frame, textvariable=self.plasmid_message_var)
        self.plasmid_message_label.pack(pady=5)

        self.input_frame = ttk.Frame(self.main_frame)
        self.input_frame.pack(pady=5)

        ttk.Label(self.input_frame, text="Plasmid Name:").grid(row=0, column=0, padx=5, pady=5)
        self.plasmid_name_entry = ttk.Entry(self.input_frame)
        self.plasmid_name_entry.grid(row=0, column=1, padx=5, pady=5)

        ttk.Label(self.input_frame, text="Plasmid ID:").grid(row=1, column=0, padx=5, pady=5)
        self.plasmid_id_entry = ttk.Entry(self.input_frame)
        self.plasmid_id_entry.grid(row=1, column=1, padx=5, pady=5)

        ttk.Label(self.input_frame, text="Plasmid Description:").grid(row=2, column=0, padx=5, pady=5)
        self.plasmid_description_entry = ttk.Entry(self.input_frame, font = ("Helvetica" , 12))
        self.plasmid_description_entry.grid(row=2, column=1, padx=5, pady=5)

        self.submit_button = ttk.Button(self.input_frame, text="Submit", command=self.submit_plasmid)
        self.submit_button.grid(row=3, column=0, columnspan=2, pady=5)

        self.input_frame.grid_columnconfigure(0, weight=1)
        self.input_frame.grid_columnconfigure(1, weight=1)

        self.input_frame.pack_forget()
        self.insertion_button.pack_forget()
        self.csv_entry.pack_forget()
        self.csv_label.pack_forget()

   

    def load_plasmids(self):
       
   
        data = self.library._retrieve_all_json_plasmids()
        for unique_id, plasmid_json in data:
            plasmid_data = json.loads(plasmid_json)
            plasmid_id = plasmid_data["unique_id"]
            plasmid_name = plasmid_data["name"]
            
            for module in plasmid_data["structure"]:
                if module["module_type"] == "cargo":
                    cargo_id = module["unique_id"]
                    cargo_name = module["name"]
                if module["module_type"] == "marker":
                    marker_id = module["unique_id"]
                    marker_name = module["name"]
                if module["module_type"] == "replication":
                    replication_id = module["unique_id"]
                    replication_name = module["name"]

            if self.tree.exists(plasmid_id):
                self.tree.item(plasmid_id, values=(plasmid_name, plasmid_id, cargo_name, cargo_id, marker_name, marker_id, replication_name, replication_id))
            else:
                self.tree.insert("", "end", iid=plasmid_id, values=(plasmid_name, plasmid_id, cargo_name, cargo_id, marker_name, marker_id, replication_name, replication_id))

    

    def select(self, event):
        self.selection_value = self.selection.get()
        print(self.selection_value)

    def insert_plasmids(self):
        csv_name = self.csv_entry.get()
        module_to_plasmid = self.library._module_to_plasmid()

        if self.module_type == "cargo":
            if self.selection_value: #Just for single/multi-channel 
                for pairs in self.plasmid_objects: #Just do add plasmid_id to the plasmid_object list for reference
                    donor_plasmid = pairs[0]
                    recipient_plasmid = pairs[1]
                    desired_plasmid_module_list = (donor_plasmid.structure[0].unique_id,recipient_plasmid.structure[2].unique_id,recipient_plasmid.structure[4].unique_id)
                    desired_plasmid_id = module_to_plasmid[desired_plasmid_module_list]
                    pairs.append(desired_plasmid_id)

        if self.module_type == "marker":
            if self.selection_value: #Just for single/multi-channel 
                for pairs in self.plasmid_objects:
                    donor_plasmid = pairs[0]
                    recipient_plasmid = pairs[1]
                    desired_plasmid_module_list = (recipient_plasmid.structure[0].unique_id,donor_plasmid.structure[2].unique_id,recipient_plasmid.structure[4].unique_id)
                    desired_plasmid_id = module_to_plasmid[desired_plasmid_module_list]
                    pairs.append(desired_plasmid_id)

        if self.module_type == "replication":
            if self.selection_value: #Just for single/multi-channel 
                for pairs in self.plasmid_objects:
                    donor_plasmid = pairs[0]
                    recipient_plasmid = pairs[1]
                    desired_plasmid_module_list = (recipient_plasmid.structure[0].unique_id,recipient_plasmid.structure[2].unique_id,donor_plasmid.structure[4].unique_id)
                    desired_plasmid_id = module_to_plasmid[desired_plasmid_module_list]
                    pairs.append(desired_plasmid_id)
        
        
        if self.selection_value == "Single Channel":
            input_csv = process_plasmid_list(self.plasmid_objects,self.module_type)  
            import_single_parameters(csv_name,input_csv)
            import_parameters_ligation(csv_name,input_csv,"plasmid")
        elif self.selection_value =="Multi-Channel":
            input_csv = process_plasmid_list(self.plasmid_objects,self.module_type)
            import_multi_parameters(csv_name,input_csv)
            import_parameters_ligation(csv_name,input_csv,"plasmid")
        elif self.selection_value == None:
            print("Please select a protocol type")




    def submit_plasmid(self):
        plasmid_name = self.plasmid_name_entry.get()
        plasmid_id = self.plasmid_id_entry.get()
        plasmid_description = self.plasmid_description_entry.get()

        if self.module_type == "cargo":
            plasmid = Plasmid(plasmid_name, plasmid_id, plasmid_description)
            plasmid.structure = [self.plasmid_objects[-1][0].structure[0],self.plasmid_objects[-1][1].structure[2],self.plasmid_objects[-1][1].structure[4]]
        
        if self.module_type == "marker":
            plasmid = Plasmid(plasmid_name, plasmid_id, plasmid_description)
            plasmid.structure = [self.plasmid_objects[-1][1].structure[0],self.plasmid_objects[-1][0].structure[2],self.plasmid_objects[-1][1].structure[4]]

        if self.module_type == "replication":
            plasmid = Plasmid(plasmid_name, plasmid_id, plasmid_description)
            plasmid.structure = [self.plasmid_objects[-1][1].structure[0],self.plasmid_objects[-1][1].structure[2],self.plasmid_objects[-1][0].structure[4]]
        
        if self.plasmid_objects[-1][1].scar_info:
            scar_data = self.plasmid_objects[-1][1].scar_info
            for scar_id, scar_info in scar_data.items():
                scar_part = Part(
                    name=scar_info['name'],
                    unique_id=scar_id,
                    sequence=scar_info['sequence'],
                    role='scar',
                    description=scar_info['description']
                )
                plasmid.insert_scar(scar_part, scar_info['module'], scar_info['position'])

        checker = self.library._constructor_add_plasmid(plasmid)

        if checker:
            self.input_frame.pack_forget()
            self.select_button.pack(pady=5)
            self.delete_button.pack(pady=5)
            self.plasmid_message_var.set("")
            self.csv_label.pack(pady=5)
            self.csv_entry.pack(pady = 5)
            self.insertion_button.pack(pady=5)
            self.message_var.set("                                Now select the next Donor Plasmid \n If you have finished your selections, click the 'Insert Plasmid' button")
            self.clear_input_fields()


    def clear_input_fields(self):
        self.plasmid_name_entry.delete(0, tk.END)
        self.plasmid_id_entry.delete(0, tk.END)
        self.plasmid_description_entry.delete(0, tk.END)


    def select_plasmid(self):
        selected_items = self.tree.selection() 
        unique_id  = str(selected_items[0])
        plasmid =  self.library.retrieve_plasmid(unique_id)
        self.temp_selection.append(plasmid)
        if len(self.temp_selection) == 2:
            if self.module_type == "cargo":
                module_id_list = [self.temp_selection[0].structure[0].unique_id,self.temp_selection[1].structure[2].unique_id,self.temp_selection[1].structure[4].unique_id]
                module_id_check = self.library._retrieve_all_module_id()
                self.csv_label.pack(pady=5)
                self.csv_entry.pack(pady = 5)
                self.insertion_button.pack(pady=5)

            if self.module_type == "marker":
            
                module_id_list = [self.temp_selection[1].structure[0].unique_id,self.temp_selection[0].structure[2].unique_id,self.temp_selection[1].structure[4].unique_id]
                module_id_check = self.library._retrieve_all_module_id()
                self.csv_label.pack(pady=5)
                self.csv_entry.pack(pady = 5)
                self.insertion_button.pack(pady=5)
        
            if self.module_type == "replication":
           
                module_id_list = [self.temp_selection[1].structure[0].unique_id,self.temp_selection[1].structure[2].unique_id,self.temp_selection[0].structure[4].unique_id]
                module_id_check = self.library._retrieve_all_module_id()
                self.csv_label.pack(pady=5)
                self.csv_entry.pack(pady = 5)
                self.insertion_button.pack(pady=5)



            if module_id_list not in module_id_check:
                self.input_frame.pack(pady=5)
                self.select_button.pack_forget()
                self.delete_button.pack_forget()
                self.plasmid_message_var.set("Generated Plasmid is new, please insert details below to add to database")
                self.insertion_button.pack_forget()
                self.csv_entry.pack_forget()
                self.csv_label.pack_forget()
                self.message_var.set("")

            self.plasmid_objects.append(self.temp_selection)
            self.message_var.set("                                Now select the next Donor Plasmid \n If you have finished your selections, click the 'Insert Plasmid' button")
            self.delete_message_var.set("")
            self.temp_selection = []
            print(f"Donor Plasmid: {self.plasmid_objects[-1][-2].name} and Recipient Plasmid: {self.plasmid_objects[-1][-1].name} have been added")

        else:
            self.message_var.set("Please select the Recipient Plasmid")
            self.delete_message_var.set("")
            self.insertion_button.pack_forget()
            self.csv_entry.pack_forget()
            self.csv_label.pack_forget()
            self.select_button.pack(pady=5)
            self.delete_button.pack(pady=5)

    def remove_last_selected(self):
        if len(self.temp_selection) == 1:
            self.delete_message_var.set(f"Donor Plasmid: {self.temp_selection[-1].name} has been removed")
            self.temp_selection.pop()
            self.message_var.set("                                Now select the next Donor Plasmid \n If you have finished your selections, click the 'Insert Plasmid' button")
            
        
        elif self.plasmid_objects:
            self.delete_message_var.set(f"Donor Plasmid: {self.plasmid_objects[-1][-2].name} and Recipient Plasmid: {self.plasmid_objects[-1][-1].name} have been removed")
            self.plasmid_objects.pop()
            self.message_var.set("                                Now select the next Donor Plasmid \n If you have finished your selections, click the 'Insert Plasmid' button")

    def print_selected(self):
        print("The Plasmids which have been selected:")
        for i in self.plasmid_objects:
            print(f"Donor Plasmid: {i[0].name}, Recipient Plasmid: {i[1].name}")



    def run(self):
        self.root.mainloop()


def plasmid_insertion_tool(library,module):  
    builder = PlasmidConstruct(library,module) 
    builder.run()

