from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio import SeqIO
import re

class Part:
    """Class represents the most fundamental componenet of a plasmid:
    Data attributes:
    Name: String
    Unique ID: String
    Sequence: String (only ACTG and no spaces)
    Role: String, gives function for further info and serialiastion
    Description: String, add some information about the part"""

    def __init__(self, name, unique_id, sequence, role, description=""):
        self.name = name
        self.unique_id = unique_id
        self.sequence = sequence
        self.role = role
        self.description = description
        self._dynamic = True

    @property
    def name(self):
        """gets the name attribute"""
        return self._name

    @name.setter
    def name(self, name):
        """sets the name attribute, type: String"""
        if not isinstance(name, str):
            raise TypeError("Name must be a string")
        self._name = name

    @property
    def unique_id(self):
        """gets the unique_id attribute"""
        return self._unique_id

    @unique_id.setter
    def unique_id(self, id):
        """gets the name attribute, type: String"""
        if not isinstance(id, str):
            raise TypeError("ID must be a string")
        self._unique_id = id

    @property
    def sequence(self):
        """gets the sequence attribute"""
        return self._sequence

    @sequence.setter
    def sequence(self, sequence):
        """sets the sequence attribute, type: String, Only contains A,C,T or G with no spaces"""
        if not isinstance(sequence, str):
            raise TypeError("Name must be a string")
        if not all(c in "ACTG" for c in sequence) or " " in sequence:
            raise ValueError("Sequence must only contain A,C,T,G with no spaces")
        self._sequence = sequence

    @property
    def role(self):
        """gets the role attribute"""
        return self._role

    @role.setter
    def role(self, role):
        """sets the role attribute, type: String, must only be certain roles also:
            "gene",
            "abr",
            "ori",
            "terminator",
            "promoter",
            "scar",
            "re",
            "orit",
            "cargo",
            "enhancer",
            "attenuator"
            "misc",
            "scar",
            "operator",
            "rbs"
            """
        

        if not isinstance(role, str):
            raise TypeError("Role must be a string")
        if role not in [
            "gene",
            "abr",
            "ori",
            "terminator",
            "promoter",
            "scar",
            "re",
            "orit",
            "cargo",
            "enhancer",
            "attenuator",
            "misc",
            "scar",
            "operator",
            "rbs",
            "regulatory"
        ]:
            raise ValueError(
                "Role must be either gene, abr, ori, terminator, promoter, scar, re, orit, cargo, enhancer, attenuator,misc,scar,operator, rbs, regulatory"
            )

        self._role = role

    @property
    def description(self):
        """gets the description attribute"""
        return self._description

    @description.setter
    def description(self, description):
        """sets the description attribute, type: String"""
        if not isinstance(description, str):
            raise TypeError("Description must be a string")
        self._description = description

    def to_dict(self):
        return {
            'name': self.name,
            'unique_id': self.unique_id,
            'sequence': self.sequence,
            'role': self.role,
            'description': self.description
        }


class Module:
    """Class represents the most module componenet of a plasmid:
    Data attributes:
    Name: String
    Unique ID: String
    Module Type: String (either cargo, marker, replication)
    Structure: List of Part objects
    Description: String, add some information about the part
    Fetched: True or False, allows for different structure construction if sequence is imported
    """


    def __init__(
        self, name, unique_id, module_type, description="", structure=[]
    ):
        self.name = name
        self.unique_id = unique_id
        self.module_type = module_type
        self.description = description
        self.structure = structure
        self._dynamic = True

    @property
    def name(self):
        """gets the name attribute"""
        return self._name

    @name.setter
    def name(self, name):
        """sets the name attribute, type: String"""
        if not isinstance(name, str):
            raise TypeError("Name must be a string")
        self._name = name

    @property
    def unique_id(self):
        """gets the unique_id attribute"""
        return self._unique_id

    @unique_id.setter
    def unique_id(self, id):
        """sets the unique_id attribute, type: String"""
        if not isinstance(id, str):
            raise TypeError("ID must be a string")
        self._unique_id = id 
        
      
        
    @property
    def module_type(self):
        """gets the module_type attribute"""
        return self._module_type

    @module_type.setter
    def module_type(self, module):
        """sets the unique_id attribute"""
        if module not in ["cargo", "marker", "replication", "invariant"]:
            raise ValueError("Module must be either a cargo, marker, replication")
        self._module_type = module


    @property
    def description(self):
        """gets the description attribute"""
        return self._description

    @description.setter
    def description(self, description):
        """sets the description attribute, type: String"""
        if not isinstance(description, str):
            raise TypeError("Description must be a string")
        self._description = description

    @property
    def structure(self):
        """gets the structure attribute"""
        return self._structure

    @structure.setter
    def structure(self, structure_list):
        """Function which sets the structure
        Input: list of part(s) objects
        Depending on module_type a predefined structure with corresponding flanking RE sites are assigned
        """

        if not isinstance(structure_list, list):
            raise TypeError("Input structure must be of Type: list")
        for parts in structure_list:
            if not isinstance(parts, Part):
                raise TypeError("List should contain Part types only")

        if (
            self.module_type == "cargo"
        ):  # Creates Part objects for flanking RE sites dynamically
            self.cargo = None
            self.PacI = Part("PacI", "paci", "TTAATTAA", "re")
            self.PacI._dynamic = False
            self.SpeI = Part("SpeI", "spei", "ACTAGT", "re")
            self.SpeI._dynamic = False
            self._structure = [self.PacI, self.cargo, self.SpeI]
            # The input part list is inserted by a slice to remove the placeholder
            self._structure[1:2] = (
                structure_list  
            )

        if self.module_type == "replication":
            self.replication = None
            self.FseI = Part("FseI", "fsei", "GGCCGGCC", "re")
            self.FseI._dynamic = False
            self.AscI = Part("AscI", "asci", "GGCGCGCC", "re")
            self.AscI._dynamic = False
            self._structure = [self.FseI, self.replication, self.AscI]
            # The input part list is inserted by a slice to remove the placeholder
            self._structure[1:2] = (
                structure_list  
            )

        if self.module_type == "marker":
            if structure_list[0].sequence != "ATTTAAAT":
                raise ValueError("Please ensure SwaI part is at the start of the structure with the correct sequence")
            
            pshai_sequence = structure_list[-1].sequence
            if ((pshai_sequence[:3] + pshai_sequence[-3:]) != "GACGTC"):
                raise ValueError("Please ensure PshAI part is at the start of the structure with the correct sequence")

            self.SwaI = structure_list[0]
            self.PshAI = structure_list[-1]
            self._structure = structure_list

       
        if self.module_type == "invariant":
            self._structure = structure_list




    def get_sequence(self):
        """Obtains a sequence by concatenating sequence of parts together"""
        sequence = ""
        for i in self.structure:
            sequence += i.sequence
        return sequence
    
    def to_dict(self):
        return {
            'name': self.name,
            'unique_id': self.unique_id,
            'module_type': self.module_type,
            'description': self.description,
            'structure': [part.to_dict() for part in self.structure if part._dynamic]
        }


class Plasmid:
    """Class which builds final plasmid:
    Data attributes:
    Name: String
    Unique ID: String
    Structure: List of Module objects
    Description: String, add some information about the part
    """

    def __init__(self, name, unique_id, description=""):
        self.name = name
        self.unique_id = unique_id
        self.description = description
        self.scar_info = {}

        # Three placeholder variables to just provide template structure
        self.cargo_module = None
        self.marker_module = None  
        self.replication_module = None

        # Creating the non-variable modules on construction of an object
        self.t1_part = Part(
            "t1_part",
            "t1_part",
            "CAGCTGTCTAGGGCGGCGGATTTGTCCTACTCAGGAGAGCGTTCACCGACAAACAACAGATAAAACGAAAGGCCCAGTCTTTCGACTGAGCCTTTCGTTTTATTTGATGCCT",
            "terminator",
        )
        self.t1 = Module("t1", "t1_module", "invariant", structure=[self.t1_part])
        self.t1._dynamic = False
        self.t1_part._dynamic = False
        self.t0_part = Part(
            "t0_part",
            "t0_part",
            "CTTGGACTCCTGTTGATAGATCCAGTAATGACCTCAGAACTCCATCTGGATTTGTTCAGAACGCTCGGTTGCCGCCGGGCGTTTTTTATTGGTGAGAATCCAG",
            "terminator",
        )
        self.SanDI = Part("SanDI_part", "sandi", "GGGTCCC", "re")
        self.t0 = Module(
            "t0", "t0_module", "invariant", structure=[self.t0_part, self.SanDI]
        )
        self.t0._dynamic = False
        self.t0_part._dynamic = False
        self.SanDI._dynamic = False
        self.orit_part = Part(
            "orit_part",
            "orit_part",
            "CTTTTCCGCTGCATAACCCTGCTTCGGGGTCATTATAGCGATTTTTTCGGTATATCCATCCTTTTTCGCACGATATACAGGATTTTGCCAAAGGGTTCGTGTAGACTTTCCTTGGTGTATCCAACGGCGTCAGCCGGGCAGGATAGGTGAAGTAGGCCCACCCGCGAGCGGGTGTTCCTTCTTCACTGTCCCTTATTCGCACCTGGCGGTGCTCAACGGGAATCCTGCTCTGCGAGGCTGGCCGTA",
            "orit",
        )
        self.orit = Module(
            "orit", "orit_module", "invariant", structure=[self.orit_part]
        )
        self.orit._dynamic = False
        self.orit_part._dynamic = False
        # Skip first validation on initilisation so ._
        self._structure = [
            self.cargo_module,
            self.t0,
            self.marker_module,
            self.orit,
            self.replication_module,
            self.t1,
        ]

    @property
    def name(self):
        """Gets the name attribute"""
        return self._name

    @name.setter
    def name(self, name):
        """Sets the name attribute, type: String"""
        if not isinstance(name, str):
            raise TypeError("Name must be a string")
        self._name = name

    @property
    def unique_id(self):
        """Sets the unique_id attribute"""
        return self._unique_id

    @unique_id.setter
    def unique_id(self, id):
        """Sets the unique_id attribute, type: String"""
        if not isinstance(id, str):
            raise TypeError("ID must be a string")
        self._unique_id = id

    @property
    def description(self):
        """gets the description attribute"""
        return self._description

    @description.setter
    def description(self, description):
        """sets the description attribute, type: String"""
        if not isinstance(description, str):
            raise TypeError("Description must be a string")
        self._description = description


    @property
    def structure(self):
        """Gets the structure attribute"""
        return self._structure

    @structure.setter
    def structure(self, structure):
        """Sets the structure attribute
        Input: A list of three individual seperate modules, type Cargo,Marker,Replication"""
        if not isinstance(structure, list):
            raise TypeError("Input structure must be of Type: list")
        for modules in structure:
            if not isinstance(modules, Module):
                raise TypeError("List should contain Module types only")

        # Creates a list containing the module types and compares it with predefined structure to see if its correct
        correct_val = ["cargo", "marker", "replication"]
        order = [module.module_type for module in structure]

        if correct_val != order:
            raise ValueError(
                "Please ensure the modules are in the correct order: [Cargo,Marker,Replication]"
            )
        
        #For use in the future construction methods
        self.components = structure 

        iterVal = iter(structure)
        for parts, element in enumerate(self._structure):
            if (
                element is None
            ):  # Finds which variables are "None" in structure AKA placeholders and then swaps them for each corresponding module
                self._structure[parts] = next(iterVal)

    def get_sequence(
        self,
    ):  # it iterates through each module in structure and then calls their get_sequence function and concatenates all together
        plasmid_sequence = ""
        for module in self.structure:
            plasmid_sequence += module.get_sequence()
        return plasmid_sequence
    

    def insert_scar(self, scar, module=None, position=None):
        """ "Method to insert scars:
        Modules: T1,T0,OriT
        Position: before or after"""
        if not isinstance(scar,Part):
            raise(TypeError("Scar must be of type:"))

        module_dict = {"T1": self.t1, "T0": self.t0, "OriT": self.orit}

        module_obj = module_dict.get(
            module
        )  # Get the module object based on the module name

        if (
            module
        ):  # Inserts the scar depending on if it comes before or after the module
            if position == "before":
                module_obj.structure.insert(0, scar)
            elif position == "after":
                module_obj.structure.append(scar)

        self.scar_info[scar.unique_id] = {
                'name': scar.name,
                'sequence': scar.sequence,
                'description': scar.description,
                'module': module,
                'position': position
            }

    def gather_data(plasmid):
        """Method which gathers the data of each part in the plasmid and puts it all into a dictionary with key being name and values being other data attributes and locations"""
        current_position = 0
        my_dict = {}

        for module in plasmid.structure:
            for part in module.structure:
                # Calculate the start position based on current position in the sequence
                start_position = current_position
                end_position = start_position + len(part.sequence) - 1

                my_dict[part.name] = [start_position, part.unique_id, part.role, part.description, end_position]

                # Update current position
                current_position += len(part.sequence)

        return my_dict
    

    def to_dict(self):
        return {
            'name': self.name,
            'unique_id': self.unique_id,
            'description': self.description,
            'structure': [module.to_dict() for module in self.structure if module._dynamic],
            'scars' : self.scar_info
        }
    



