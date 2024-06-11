from opentrons import protocol_api
import csv
import io 
# requirements
requirements = {"robotType": "OT-2", "apiLevel": "2.16"}

input_notepad_paramaters_here = """
PASTE CSV HERE AND DELETE THIS MESSAGE
"""


def get_paramaters(csv_data):
    normalised_data = "\n".join([line.strip() for line in csv_data.strip().splitlines() if line.strip()])

    data = []
    csvfile = io.StringIO(normalised_data)
    reader = csv.reader(csvfile)
    plasmid_count = 0
    data = []
    
    for row in reader:
    # Check if it's a "Plasmid Name" row
        if "Plasmid Name" in row[0]:
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
                data.append([rows[(6*i)+1][q],rows[(6*i)+2][q],rows[(6*i)+5][q]])

        for q in range((full_cols*2) +1,(full_cols*2) +3):
            for i in range(0,int(remaining_cells/2)):
                data.append([rows[(6*i)+1][q],rows[(6*i)+2][q],rows[(6*i)+5][q]])
    if full_cols == 0:
        for q in range(1,3):
            for i in range(0,int(remaining_cells/2)):
                data.append([rows[(6*i)+1][q],rows[(6*i)+2][q],rows[(6*i)+5][q]])


    data_list = data
    
    def sort_custom_simple(data):
    # Simple sorting function
        sorted_data = sorted(data, key=lambda x: (x[0][0], int(x[0][1:])))
        return sorted_data
 
    sorted_data = sort_custom_simple(data_list)

    def group_in_pairs(data):
        return [data[i:i+2] for i in range(0, len(data), 2)]

    paired_data = group_in_pairs(sorted_data)

    return paired_data
 


paired_data = get_paramaters(input_notepad_paramaters_here)


water_vol = 0
for i in paired_data:
    water_vol += 17 - (float(i[0][2]) + float(i[1][2]))


def run(protocol: protocol_api.ProtocolContext):

    plate_sample = protocol.load_labware("opentrons_96_wellplate_200ul_pcr_full_skirt", 1)
    plate_reagents = protocol.load_labware("opentrons_96_wellplate_200ul_pcr_full_skirt", 2)
    right_tips1 = protocol.load_labware("opentrons_96_tiprack_20ul", 3)
    right_tips2 = protocol.load_labware("opentrons_96_tiprack_20ul", 4)
    tube_rack = protocol.load_labware("opentrons_24_tuberack_eppendorf_2ml_safelock_snapcap", 5)

    right_pipette = protocol.load_instrument("p20_single_gen2", "right", tip_racks=[right_tips1,right_tips2])

    t4_ligase_buffer = protocol.define_liquid(name="T4 DNA Ligase Buffer (10X)", description= "Buffer for T4 Ligase enzyme", display_color= "#FF7820")
    t4_ligase = protocol.define_liquid(name="T4 DNA Ligase", description= "T4 DNA Ligase ", display_color= "#84FBB1")
    water = protocol.define_liquid(name="Water", description= "Water", display_color= "#D7FC26")


    plate_reagents["A1"].load_liquid(liquid=t4_ligase_buffer, volume = (len(paired_data))*2)
    plate_reagents["B1"].load_liquid(liquid=t4_ligase, volume = (len(paired_data)))
    tube_rack["A1"].load_liquid(liquid=water, volume = round(water_vol) + 100) 
    
    for plasmid_pairs in paired_data:
        donor_well = plasmid_pairs[0][0]
        donor_volume = float(plasmid_pairs[0][2])
        recipient_well = plasmid_pairs[1][0]
        recipient_volume = float(plasmid_pairs[1][2])
        water_volume = 17 - donor_volume - recipient_volume 

        right_pipette.pick_up_tip()
        right_pipette.aspirate(water_volume,tube_rack["A1"])
        right_pipette.dispense(water_volume,plate_sample[donor_well])
        right_pipette.drop_tip()
        right_pipette.pick_up_tip()
        right_pipette.aspirate(2,plate_reagents["A1"])  #Pick buffer up
        right_pipette.air_gap(volume = 1)
        right_pipette.aspirate(1,plate_reagents["B1"])  #Pick ligase up
        right_pipette.air_gap(volume = 1)
        right_pipette.aspirate(recipient_volume,plate_sample[recipient_well])
        right_pipette.air_gap(volume=1)
        right_pipette.dispense(recipient_volume + 6, plate_sample[donor_well])
        right_pipette.mix(3, 0.7*20, plate_sample[donor_well])
        right_pipette.drop_tip()




