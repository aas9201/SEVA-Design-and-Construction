from plasmid import Part, Module, Plasmid
import tkinter as tk
from tkinter import ttk
import json
import textwrap
from library import Library


class ModuleMaker:
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
        self.selection_message_var = tk.StringVar()
        self.part_selection = []
        self.part_objects = []
        self.insert_input = {}
        self.part_sequences = {}
        self.selection_value = None
        self.frame = ttk.Frame(self.root)
        self.frame.grid(row=0, column=0, sticky='nsew')
        self.root.title("ModuleMaker Tool")
        self.create_widgets_frame()
        self.load_parts()

    def create_widgets_frame(self):
        instructions = ttk.Label(self.frame, text="Select and add parts to create sequence.")
        instructions.pack(pady=10)

        self.tree2 = ttk.Treeview(self.frame)
        self.tree2["columns"] = ("part_name", "part_id", "part_role", "part_description")
        self.tree2.column("#0", width=0, stretch=tk.NO)
        self.tree2.column("part_name", anchor=tk.W, width=100)
        self.tree2.column("part_id", anchor=tk.W, width=120)
        self.tree2.column("part_role", anchor=tk.W, width=150)
        self.tree2.column("part_description", anchor=tk.W, width=250)
        self.tree2.heading("#0", text="", anchor=tk.W)
        self.tree2.heading("part_name", text="Part Name", anchor=tk.W)
        self.tree2.heading("part_id", text="Part ID", anchor=tk.W)
        self.tree2.heading("part_role", text="Role", anchor=tk.W)
        self.tree2.heading("part_description", text="Description", anchor=tk.W)
        self.tree2.pack(fill=tk.BOTH, expand=True)

        add_part_button = ttk.Button(self.frame, text="Add Part", command=self.add_part)
        add_part_button.pack(pady=10)

        add_print_button = ttk.Button(self.frame, text="Print Part Sequence", command=self.print_part)
        add_print_button.pack(pady=10)


        self.delete_button = ttk.Button(self.frame, text="Remove Last Selection", command=self.remove_last_selected)
        self.delete_button.pack(pady=10)

        self.selection_message_label = ttk.Label(self.frame, textvariable=self.selection_message_var)
        self.selection_message_label.pack(pady=5)
        self.selection_message_var.set("Current Selection:")

        ttk.Label(self.frame, text="Name:").pack(padx=5, pady=5)
        self.name_entry = ttk.Entry(self.frame)
        self.name_entry.pack(pady=5)

        ttk.Label(self.frame, text="ID:").pack(padx=5, pady=5)
        self.id_entry = ttk.Entry(self.frame)
        self.id_entry.pack(pady=5)

        ttk.Label(self.frame, text="Description:").pack(padx=5, pady=5)
        self.description_entry = tk.Text(self.frame, height=3, width=25, font=('Helvetica', 10), wrap = tk.WORD)
        self.description_entry.pack(pady=5)


        options = ["Cargo", "Marker", "Replication", "MCS Insert"]
        self.selection = ttk.Combobox(self.frame, values=options, state= "readonly")
        self.selection.set("Select Construct Type")
        self.selection.bind("<<ComboboxSelected>>", self.select)
        self.selection.pack(pady=10)


        self.insert_button = ttk.Button(self.frame, text="Insert Parts", command=lambda: self.insert_parts())
        self.insert_button.pack(pady=10)
        
        


    
    def load_parts(self):
        data = self.library._retrieve_all_json_parts()
        part_id = part_name = part_role = part_description = part_sequence = ""

        for unique_id, part_json in data:
            part_data = json.loads(part_json)
            part_id = part_data["unique_id"]
            part_name = part_data["name"]
            part_role = part_data["role"]
            part_description = part_data["description"]
            self.part_sequences[unique_id] = part_data["sequence"]
            
            if self.tree2.exists(unique_id):
                self.tree2.item(unique_id, values=(part_name, part_id, part_role, part_description))
            else:
                self.tree2.insert("", "end", iid=unique_id, values=(part_name, part_id, part_role, part_description))

    

    def add_part(self):

        selected_items = self.tree2.selection() 
        unique_id  = str(selected_items[0])
        self.part_selection.append(unique_id)

        self.selection_message_var.set(f"Current Selection:{self.part_selection}")

    def print_part(self):
        selected_items = self.tree2.selection()
        if selected_items:  # Ensure there's a selected item
            unique_id = selected_items[0]  # Get the first selected item's ID
            if unique_id in self.part_sequences:
                sequence = self.part_sequences[unique_id]
                print(f"Sequence for {unique_id}: {sequence}")
            else:
                print("Sequence not found.")
        else:
            print("No part selected.") 


    def insert_parts(self):
        if self.selection_value:
            for i in self.part_selection:
                part_object = self.library.retrieve_part(i)
                self.part_objects.append(part_object)

        name = self.name_entry.get()
        id = self.id_entry.get()
        description = self.description_entry.get()

        if self.selection_value == "MCS Insert":
            insert_data = {'name': name,
            'unique_id': id,
            'description': description,
            'structure': self.part_selection}     
              
            checker = self.library._constructor_add_insert(insert_data)

            if checker:
                self.clear_input_fields()
                self.part_selection = []
                self.part_objects = []
                self.selection_message_var.set(f"Current Selection:{self.part_selection}")
        
        elif self.selection_value == "Cargo":
            cargo_module = Module(name,id,"cargo",description, structure=self.part_objects)
            checker = self.library._constructor_add_module(cargo_module)
            if checker:
                self.clear_input_fields()
                self.part_selection = []
                self.part_objects = []
                self.selection_message_var.set(f"Current Selection:{self.part_selection}")

        elif self.selection_value == "Marker":
            marker_module = Module(name,id,"marker",description, structure=self.part_objects)
            checker = self.library._constructor_add_module(marker_module)
            if checker:
                self.clear_input_fields()
                self.part_selection = []
                self.part_objects = []
                self.selection_message_var.set(f"Current Selection:{self.part_selection}")

        elif self.selection_value == "Replication":
            replication_module = Module(name,id,"replication",description, structure=self.part_objects)
            checker = self.library._constructor_add_module(replication_module)
            if checker:
                self.clear_input_fields()
                self.part_selection = []
                self.part_objects = []
                self.selection_message_var.set(f"Current Selection:{self.part_selection}")



    def select(self, event):
        self.selection_value = self.selection.get()
        print(self.selection_value)

    def clear_input_fields(self):
        self.name_entry.delete(0, tk.END)
        self.id_entry.delete(0, tk.END)
        self.description_entry.delete(0, tk.END)


    def remove_last_selected(self):
        if self.part_selection:
            self.part_selection.pop()
            self.selection_message_var.set(f"Current Selection:{self.part_selection}")
            
            

    def run(self):
        self.root.mainloop()

def module_insertion_tool(library):  
    builder = ModuleMaker(library) 
    builder.run()  


tester1= Library("example_database")
module_insertion_tool(tester1)