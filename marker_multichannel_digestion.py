from opentrons import protocol_api
import math
import csv
import io
# requirements
requirements = {"robotType": "OT-2", "apiLevel": "2.16"}

metadata = {
    "protocolName": "Multichannel Restriction Enzyme Digestion",
    "author": "Adnan Salem",
    "description": "Automated restriction enzyme digestion protocol for initial plasmid cloning steps, utilising dual pipetting systems for precision handling of enzymes, buffer, and water across multiple samples. Optimised for consistent preparation in 96-well plate format ",
}

import_csv_notepad_data_here ="""
PASTE CSV HERE AND DELETE THIS MESSAGE
"""


def get_multi_parameter(csv_data):

    normalised_data = "\n".join([line.strip() for line in csv_data.strip().splitlines() if line.strip()])
    csvfile = io.StringIO(normalised_data)
    data = []
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
    

    csvfile.seek(0)
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
    
    return data_list



param = get_multi_parameter(import_csv_notepad_data_here)

plasmid_no = len(param) - 1
column_no = math.ceil(plasmid_no/8)

print(column_no)
def run(protocol: protocol_api.ProtocolContext):

    enzyme_volume1 = int(param[0][0])
    enzyme_volume2 = int(param[0][1])
    plasmid_volume = int(param[0][2])
    reaction_volume1 = int(param[0][3])
    reaction_volume2 = int(param[0][4])
    plasmid_no = len(param) - 1
    column_no = math.ceil(plasmid_no/8)


    temp_mod = protocol.load_module("temperature module gen2", 1)
    temp_mod_adapter = temp_mod.load_adapter("opentrons_96_well_aluminum_block")
    temp_plate =temp_mod_adapter.load_labware("opentrons_96_wellplate_200ul_pcr_full_skirt")
    reagent_plate = protocol.load_labware("opentrons_96_wellplate_200ul_pcr_full_skirt", 2)
    reservoir = protocol.load_labware("nest_12_reservoir_15ml", 5)
    right_tips1 = protocol.load_labware("opentrons_96_tiprack_20ul", 7)
    right_tips2 = protocol.load_labware("opentrons_96_tiprack_20ul", 9)
    left_tips1 = protocol.load_labware("opentrons_96_tiprack_300ul", 8)
    left_tips2 = protocol.load_labware("opentrons_96_tiprack_300ul", 11)

    left_pipette = protocol.load_instrument("p300_multi_gen2", "left", tip_racks=[left_tips1,left_tips2])
    right_pipette = protocol.load_instrument("p20_multi_gen2", "right", tip_racks=[right_tips1,right_tips2])

    
    PshAI = protocol.define_liquid(name="PshAI enzyme", description= "Enzyme for digesion", display_color= "#FF7820")
    SwaI = protocol.define_liquid(name="SwaI enzyme", description= "Enzyme for digesion", display_color= "#84FBB1")
    Buffer = protocol.define_liquid(name="Buffer", description= "Buffer", display_color= "#D7FC26")
    Water = protocol.define_liquid(name="Water", description= "Water", display_color= "#51F360")
    NaCl = protocol.define_liquid(name="NaCl", description= "1M NaCl Solution", display_color= "#007BFF")
    reagent_columns = reagent_plate.rows()[0]
    sample_columns = temp_plate.rows()[0]
    
    


    if plasmid_no % 16 == 0:
        water_multiplier = plasmid_no
    else:
        water_multiplier = (plasmid_no + (16 - plasmid_no % 16))/2

    water_no = reaction_volume1 + reaction_volume2 - 2*enzyme_volume1 - 2*enzyme_volume2 - 2*plasmid_volume - (reaction_volume1+reaction_volume2) / 10
    nacl_volume = ((reaction_volume1 + reaction_volume2) + (reaction_volume1+reaction_volume2)/10) / 10
    for i in range(0,8):
        reagent_plate.columns()[0][i].load_liquid(liquid=PshAI, volume = ((enzyme_volume1 + enzyme_volume2)*column_no)/8)
        reagent_plate.columns()[1][i].load_liquid(liquid=SwaI, volume = ((enzyme_volume1 + enzyme_volume2)*column_no)/8)
        reagent_plate.columns()[2][i].load_liquid(liquid=Buffer, volume = round(((reaction_volume1+reaction_volume2)/10)*(column_no/2)/8))
    reservoir['A1'].load_liquid(liquid=Water, volume=water_no * water_multiplier)
    reservoir['A2'].load_liquid(liquid=NaCl, volume=nacl_volume * water_multiplier)

    print(reagent_columns)
    temp_mod.set_temperature(celsius=37)
    for i in range(0,column_no):
        enzyme_volume = enzyme_volume1 if i % 2 == 0 else enzyme_volume2
        reaction_volume = reaction_volume1 if i % 2 == 0 else reaction_volume2
        buffer_volume = reaction_volume / 10
        water_volume = reaction_volume - plasmid_volume - 2*enzyme_volume - buffer_volume

        left_pipette.pick_up_tip()
        right_pipette.pick_up_tip()
        left_pipette.aspirate(water_volume, reservoir["A1"])
        left_pipette.air_gap(volume=10)
        right_pipette.aspirate(buffer_volume, reagent_columns[2])
        right_pipette.air_gap(volume=2.5)
        right_pipette.aspirate(enzyme_volume, reagent_columns[0])
        right_pipette.air_gap(volume=2.5)
        right_pipette.touch_tip()
        left_pipette.dispense(water_volume + 10, sample_columns[i])
        left_pipette.touch_tip()
        right_pipette.dispense((buffer_volume + enzyme_volume + 5), sample_columns[i])
        right_pipette.drop_tip()
        left_pipette.drop_tip()
    protocol.delay(minutes = 15)
    temp_mod.set_temperature(celsius=25)

    for i in range(0,column_no):
        enzyme_volume = enzyme_volume1 if i % 2 == 0 else enzyme_volume2
        reaction_volume = reaction_volume1 if i % 2 == 0 else reaction_volume2
        buffer_volume = reaction_volume / 10
        water_volume = reaction_volume - plasmid_volume - 2*enzyme_volume - buffer_volume
        nacl_volume = (reaction_volume + (reaction_volume/10)) / 10

        left_pipette.pick_up_tip()
        right_pipette.pick_up_tip()
        left_pipette.aspirate(nacl_volume, reservoir["A2"])
        left_pipette.dispense(nacl_volume, sample_columns[i])
        right_pipette.aspirate(enzyme_volume, reagent_columns[1])
        right_pipette.touch_tip()
        right_pipette.dispense(enzyme_volume, sample_columns[i])
        left_pipette.mix(3, 0.7*reaction_volume, sample_columns[i])
        left_pipette.drop_tip()
        right_pipette.drop_tip()
    temp_mod.deactivate()