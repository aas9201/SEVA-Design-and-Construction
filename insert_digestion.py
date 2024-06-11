from opentrons import protocol_api
import csv
import io 
# requirements
requirements = {"robotType": "OT-2", "apiLevel": "2.16"}

input_notepad_paramaters_here = """
PASTE CSV HERE AND DELETE THIS MESSAGE
"""







def get_single_parameters(csv_data):
    normalised_data = "\n".join([line.strip() for line in csv_data.strip().splitlines() if line.strip()])

    data = []
    csvfile = io.StringIO(normalised_data)
    reader = csv.reader(csvfile)
    plasmid_count = 0
    data = []
    
    for row in reader:
    # Check if it's a "Plasmid Name" row
        if "Name" in row[0]:
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
 

    if full_cols > 0:
        for q in range(1,(full_cols*2)+1):
            for i in range(0,8):
                data.append([rows[(10*i)+1][q],rows[(10*i)+2][q],rows[(10*i)+5][q],rows[(10*i)+6][q],rows[(10*i)+7][q],rows[(10*i)+8][q],rows[(10*i)+9][q]])

        for q in range((full_cols*2) +1,(full_cols*2) +3):
            for i in range(0,int(remaining_cells/2)):
                data.append([rows[(10*i)+1][q],rows[(10*i)+2][q],rows[(10*i)+5][q],rows[(10*i)+6][q],rows[(10*i)+7][q],rows[(10*i)+8][q],rows[(10*i)+9][q]])
    if full_cols == 0:
        for q in range(1,3):
            for i in range(0,int(remaining_cells/2)):
                data.append([rows[(10*i)+1][q],rows[(10*i)+2][q],rows[(10*i)+5][q],rows[(10*i)+6][q],rows[(10*i)+7][q],rows[(10*i)+8][q],rows[(10*i)+9][q]])


    data_list = data
    
    return data_list

parameters = get_single_parameters(input_notepad_paramaters_here)
print(parameters)

enzyme_counter = {}

for entry in parameters:
    enzyme_volume = int(entry[3])
    enzymes = entry[5:]  

    for enzyme in enzymes:
        if enzyme in enzyme_counter:
            enzyme_counter[enzyme] += enzyme_volume
        else:
            enzyme_counter[enzyme] = enzyme_volume
for i in enzyme_counter:
    print(enzyme_counter[i])

total_water_volume = 0
total_buffer_volume = 0

for well, name, plasmid_volume ,enzyme_volume, reaction_volume, enzyme1, enzyme2 in parameters:
    buffer_volume = int(reaction_volume) / 10
    water_volume = int(reaction_volume) - int(plasmid_volume) - 2*int(enzyme_volume) - buffer_volume
    total_water_volume += water_volume
    total_buffer_volume +=buffer_volume

def run(protocol: protocol_api.ProtocolContext):


    temp_mod = protocol.load_module("temperature module gen2", 1)
    temp_mod_adapter = temp_mod.load_adapter("opentrons_96_well_aluminum_block")
    temp_plate =temp_mod_adapter.load_labware("opentrons_96_wellplate_200ul_pcr_full_skirt")
    plate_reagents = protocol.load_labware("opentrons_96_wellplate_200ul_pcr_full_skirt", 2)
    right_tips1 = protocol.load_labware("opentrons_96_tiprack_20ul", 3)
    right_tips2 = protocol.load_labware("opentrons_96_tiprack_20ul", 4)
    reservoir = protocol.load_labware("nest_12_reservoir_15ml", 5)
    left_tips1 = protocol.load_labware("opentrons_96_tiprack_300ul", 6)
    left_tips2 = protocol.load_labware("opentrons_96_tiprack_300ul", 7)

    left_pipette = protocol.load_instrument("p300_single_gen2", "left", tip_racks=[left_tips1,left_tips2])
    right_pipette = protocol.load_instrument("p20_single_gen2", "right", tip_racks=[right_tips1,right_tips2])


    BamHI = protocol.define_liquid(name="BamHI enzyme", description="Enzyme for digestion", display_color="#D62728")
    EcoRI = protocol.define_liquid(name="EcoRI enzyme", description="Enzyme for digestion", display_color="#1E90FF")
    SacI = protocol.define_liquid(name="SacI enzyme", description="Enzyme for digestion", display_color="#228B22")
    SalI = protocol.define_liquid(name="SalI enzyme", description="Enzyme for digestion", display_color="#FF7F0E")
    HindIII = protocol.define_liquid(name="HindIII enzyme", description="Enzyme for digestion", display_color="#ECD540")
    PstI = protocol.define_liquid(name="PstI enzyme", description="Enzyme for digestion", display_color="#9467BD")
    Buffer = protocol.define_liquid(name="Buffer", description="Buffer", display_color="#E377C2")
    Water = protocol.define_liquid(name="Water", description="Water", display_color="#7F7F7F")



    plate_reagents["A1"].load_liquid(liquid=BamHI, volume = enzyme_counter["BamHI"])
    plate_reagents["A2"].load_liquid(liquid=EcoRI, volume = enzyme_counter["EcoRI"])
    plate_reagents["A3"].load_liquid(liquid=SacI, volume = enzyme_counter["SacI"])
    plate_reagents["B1"].load_liquid(liquid=SalI, volume = enzyme_counter["SalI"])
    plate_reagents["B2"].load_liquid(liquid=HindIII, volume = enzyme_counter["HindIII"])
    plate_reagents["B3"].load_liquid(liquid=PstI, volume = enzyme_counter["PstI"])
    plate_reagents["C1"].load_liquid(liquid=Buffer, volume = total_buffer_volume)
    reservoir["A1"].load_liquid(liquid=Water, volume = total_water_volume)


    temp_mod.set_temperature(celsius=37)
    for well, name, plasmid_volume ,enzyme_volume, reaction_volume, enzyme1, enzyme2 in parameters:
        plasmid_volume = int(plasmid_volume)
        enzyme_volume = int(enzyme_volume)
        reaction_volume = int(reaction_volume)
        buffer_volume = reaction_volume / 10
        water_volume = reaction_volume - plasmid_volume - 2*enzyme_volume - buffer_volume

        enzyme_to_well = {
        "BamHI": "A1",
        "EcoRI": "A2",
        "SacI": "A3",
        "SalI": "B1",
        "HindIII": "B2",
        "PstI": "B3"
        }

    
        #Donor Plasmid Row
        left_pipette.pick_up_tip()
        right_pipette.pick_up_tip()
        left_pipette.aspirate(water_volume, reservoir["A1"])
        left_pipette.air_gap(volume=10)
        right_pipette.aspirate(buffer_volume, plate_reagents["C1"])
        right_pipette.air_gap(volume=2.5)
        right_pipette.aspirate(enzyme_volume, plate_reagents[enzyme_to_well[enzyme1]])
        right_pipette.air_gap(volume=2.5)
        right_pipette.touch_tip()
        left_pipette.dispense(water_volume + 10, temp_plate[well])
        left_pipette.touch_tip()
        right_pipette.dispense((buffer_volume + enzyme_volume + 5), temp_plate[well])
        right_pipette.drop_tip()

        right_pipette.pick_up_tip()
        right_pipette.aspirate(enzyme_volume, plate_reagents[enzyme_to_well[enzyme2]])
        right_pipette.touch_tip()
        right_pipette.dispense(enzyme_volume, temp_plate[well])
        left_pipette.mix(3, 0.7*reaction_volume, temp_plate[well])
        left_pipette.drop_tip()
        right_pipette.drop_tip()

    temp_mod.deactivate()
