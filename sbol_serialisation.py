from sbol2 import ComponentDefinition, Sequence, Document, setHomespace
from sbol2 import Document
from plasmid import Plasmid, Module, Part
import os 

setHomespace("http://seva_plasmmids.com")
doc = Document()
doc.clear()


def adding_component(displayid=str, ontology=str):
    """Creates SBOL component"""
    comp = ComponentDefinition(displayid)
    comp.roles = f"http://identifiers.org/so/SO:{ontology}"
    doc.addComponentDefinition(comp)
    return comp


def assignsequence(id, code=str):
    """Assigns the sequence to SBOL Component"""
    seq = Sequence(id, code)
    return seq


def serialise_sbol(plasmid, folder_name):
    """Function which serialises in SBOL format
    Input: Plasmid object and file path"""

    current_dict = os.getcwd()
    folder_path = os.path.join(current_dict,"storage", folder_name)
    if not os.path.exists(folder_path):
            os.makedirs(folder_path)
    
    if not isinstance(plasmid, Plasmid):
        raise TypeError("Please input a plasmid object")

    doc.clear()
    primary_structure = []
    data_dict = plasmid.gather_data()
    # Dictionary to convert roles into SBOL ontology codes
    conversion = {
        "re": "0001687",  
        "cargo": "0000704",
        "abr": "0000001",
        "terminator": "0000141",
        "orit": "0000724",
        "ori": "0000296",
        "scar": "0001953",
        "promoter":"0000167",
        "gene":"0000704",
        "rbs":"0000139"
    }

    # Apply conversion
    for part, attribute in data_dict.items():
        if attribute[2] in conversion:
            attribute[2] = conversion[attribute[2]]
    # Add sequence to dictionary as need to add component and sequence in one loop for it to work
    for module in plasmid.structure:
        for part in module.structure:
            data_dict[part.name].append(part.sequence)

    for part, attibute in data_dict.items():
        component = adding_component(attibute[1], attibute[2])  # Add Sbol component
        component.sequence = assignsequence(
            component.displayId, attibute[5]
        )  # Add Sequence to component
        primary_structure.append(component)  # Add the component to primary structure

    instance_built = adding_component(
        plasmid.name, "0000155"
    )  # Create final plasmid component for hierarchy
    instance_built.assemblePrimaryStructure(
        primary_structure
    )  # Assign primary structure to plasmid component
    instance_built.compile()
    #Writes the document, full file path is needed as input for now 
    doc.write(f"{folder_path}\{plasmid.unique_id}.xml") 


