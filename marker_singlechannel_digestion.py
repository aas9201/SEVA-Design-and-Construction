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
        if "Plasmid Name" in row[0]:   #change to reaction fragment
            # Count all non-empty entries except for the first one (which is the text "Plasmid Name")
            for item in row[1:]:
                if item:  # This checks if the cell is not empty
                    plasmid_count += 1
        # print(plasmid_count)

    full_cols, remaining_cells = divmod(plasmid_count, 16)


   
    csvfile.seek(0)
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
    
    return data_list

parameters = get_single_parameters(input_notepad_paramaters_here)

for well, name, plasmid_volume ,enzyme_volume, reaction_volume in parameters:
    if well[1] == "1" or "3" or "5" or "7" or "9" or "11":
        print(well[1])

# protocol run function
def run(protocol: protocol_api.ProtocolContext):

    temp_mod = protocol.load_module("temperature module gen2", 1)
    temp_mod_adapter = temp_mod.load_adapter("opentrons_96_well_aluminum_block")
    temp_plate = temp_mod_adapter.load_labware("opentrons_96_wellplate_200ul_pcr_full_skirt")
    plate_reagents = protocol.load_labware("opentrons_96_wellplate_200ul_pcr_full_skirt", 2)
    right_tips1 = protocol.load_labware("opentrons_96_tiprack_20ul", 3)
    right_tips2 = protocol.load_labware("opentrons_96_tiprack_20ul", 4)
    tube_rack = protocol.load_labware("opentrons_24_tuberack_eppendorf_2ml_safelock_snapcap", 5)
    left_tips1 = protocol.load_labware("opentrons_96_tiprack_300ul", 6)
    left_tips2 = protocol.load_labware("opentrons_96_tiprack_300ul", 7)

    left_pipette = protocol.load_instrument("p300_single_gen2", "left", tip_racks=[left_tips1,left_tips2])
    right_pipette = protocol.load_instrument("p20_single_gen2", "right", tip_racks=[right_tips1,right_tips2])


    temp_mod.set_temperature(celsius=37)
    for well, name, plasmid_volume ,enzyme_volume, reaction_volume in parameters:
        plasmid_volume = int(plasmid_volume)
        enzyme_volume = int(enzyme_volume)
        reaction_volume = int(reaction_volume)
        buffer_volume = reaction_volume / 10
        water_volume = reaction_volume - plasmid_volume - 2*enzyme_volume - buffer_volume

        #Donor Plasmid Row
        left_pipette.pick_up_tip()
        right_pipette.pick_up_tip()
        left_pipette.aspirate(water_volume, tube_rack["A1"])
        left_pipette.air_gap(volume=10)
        right_pipette.aspirate(buffer_volume, plate_reagents["C1"])
        right_pipette.air_gap(volume=2.5)
        right_pipette.aspirate(enzyme_volume, plate_reagents["A1"])
        right_pipette.air_gap(volume=2.5)
        right_pipette.touch_tip()
        left_pipette.dispense(water_volume + 10, temp_plate[well])
        left_pipette.touch_tip()
        right_pipette.dispense((buffer_volume + enzyme_volume + 5), temp_plate[well])
        right_pipette.drop_tip()
        left_pipette.drop_tip()
    protocol.delay(minutes = 15)

    temp_mod.set_temperature(celsius=25)

    for well, name, plasmid_volume ,enzyme_volume, reaction_volume in parameters:
            plasmid_volume = int(plasmid_volume)
            enzyme_volume = int(enzyme_volume)
            reaction_volume = int(reaction_volume)
            buffer_volume = reaction_volume / 10
            water_volume = reaction_volume - plasmid_volume - 2*enzyme_volume - buffer_volume
            nacl_volume = (reaction_volume + (reaction_volume/10)) / 10

            left_pipette.pick_up_tip()
            right_pipette.pick_up_tip()
            right_pipette.aspirate(nacl_volume,tube_rack["B1"])
            right_pipette.air_gap(volume=2.5)
            right_pipette.aspirate(enzyme_volume, plate_reagents["B1"])
            right_pipette.touch_tip()
            right_pipette.dispense(enzyme_volume + nacl_volume + 2.5, temp_plate[well])
            left_pipette.mix(3, 0.7*reaction_volume, temp_plate[well])
            left_pipette.drop_tip()
            right_pipette.drop_tip()

    temp_mod.deactivate()