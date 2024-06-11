from opentrons import protocol_api

# requirements
requirements = {"robotType": "OT-2", "apiLevel": "2.16"}

# protocol run function
def run(protocol: protocol_api.ProtocolContext):

    mag_deck = protocol.load_module("magnetic module gen2", 9)
    mag_plate = mag_deck.load_labware("nest_96_wellplate_2ml_deep")

    hs = protocol.load_module("heaterShakerModuleV1", 1)
    hs_adapter = hs.load_adapter("opentrons_96_deep_well_adapter")
    hs_plate = hs_adapter.load_labware("nest_96_wellplate_2ml_deep")

    reagent_plate = protocol.load_labware("nest_96_wellplate_2ml_deep", 8)
   

    right_tips1 = protocol.load_labware("opentrons_96_tiprack_1000ul", 10)
    left_tips1 = protocol.load_labware("opentrons_96_tiprack_300ul", 11)
    left_tips2 = protocol.load_labware("opentrons_96_tiprack_300ul", 7)

    left_pipette = protocol.load_instrument("p300_multi_gen2", "left", tip_racks=[left_tips1,left_tips2])
    #need right pipette for setup  MAYBE??? Might not

    hs_column = hs_plate.rows()[0]
    mag_column = mag_plate.rows()[0]
    reagent_column = reagent_plate.rows()[0]

    
    hs.close_labware_latch()
    hs.set_target_temperature(50)
    protocol.delay(minutes = 2)
    left_pipette.transfer(300,reagent_column[8],hs_column[0])
    protocol.delay(minutes = 1) 
    left_pipette.transfer(300,reagent_column[8],hs_column[1])#1 min 
    protocol.delay(minutes=2)#1:20
    hs.set_and_wait_for_shake_speed(500) #3:20
    protocol.delay(minutes=0.5) #3:25
    hs.deactivate_shaker() 
    protocol.delay(minutes=2) #3:55
    hs.set_and_wait_for_shake_speed(500) #5:55
    protocol.delay(minutes=0.5) #6:00
    hs.deactivate_shaker() 
    protocol.delay(minutes=2) #6:30
    hs.set_and_wait_for_shake_speed(500) #8:30
    protocol.delay(minutes=0.5) #8:35
    hs.deactivate_shaker() #9:05

    left_pipette.pick_up_tip()
    left_pipette.aspirate(300,hs_column[0], rate = 0.5)
    left_pipette.dispense(300,mag_column[0], push_out=5)
    left_pipette.drop_tip()

    left_pipette.pick_up_tip()
    left_pipette.aspirate(240,reagent_column[7])
    left_pipette.air_gap(10)
    left_pipette.aspirate(50,reagent_column[6])
    left_pipette.dispense(300,mag_column[0])   
    left_pipette.mix(3,300,rate = 1.5)  #added everything into col1, 5 min wait from now
    left_pipette.drop_tip()
    hs.deactivate_heater()


    #letting col2 cool down
    left_pipette.pick_up_tip()
    left_pipette.aspirate(300,hs_column[1], rate = 0.5)
    left_pipette.dispense(300,reagent_column[11], push_out=5)
    left_pipette.drop_tip()

    left_pipette.pick_up_tip()
    left_pipette.aspirate(240,reagent_column[7])
    left_pipette.air_gap(10)
    left_pipette.aspirate(50,reagent_column[6])
    left_pipette.dispense(300,reagent_column[11])   #added everything into col2, +1 min from col1
    left_pipette.mix(3,300,rate = 1.5)
    left_pipette.drop_tip()

    protocol.delay(minutes = 4)
    mag_deck.engage()
    protocol.delay(seconds = 50)

    left_pipette.pick_up_tip()
    left_pipette.aspirate(300,reagent_column[11])
    left_pipette.dispense(300,mag_column[1])
    left_pipette.aspirate(300,reagent_column[11])
    left_pipette.dispense(300,mag_column[1])
    left_pipette.drop_tip()   #From now Col1 been on magnet for 2 min

    #First wash step of Col1
    left_pipette.pick_up_tip()
    left_pipette.aspirate(300,mag_column[0])
    left_pipette.dispense(300,reagent_column[5])     
    left_pipette.aspirate(300,mag_column[0])
    left_pipette.dispense(300,reagent_column[5])  #bin waste
    left_pipette.drop_tip()
    left_pipette.pick_up_tip()
    left_pipette.aspirate(300,reagent_column[10])  #adding wash buffer to col 1
    left_pipette.dispense(300,mag_column[0])     
    left_pipette.aspirate(300,reagent_column[10])  
    left_pipette.dispense(300,mag_column[0])
    left_pipette.drop_tip()

    #First wash step of col2
    left_pipette.pick_up_tip()
    left_pipette.aspirate(300,mag_column[1])
    left_pipette.dispense(300,reagent_column[4])     
    left_pipette.aspirate(300,mag_column[0])
    left_pipette.dispense(300,reagent_column[4])  #bin waste
    left_pipette.drop_tip()
    left_pipette.pick_up_tip()
    left_pipette.aspirate(300,reagent_column[9])  #adding wash buffer to col 2
    left_pipette.dispense(300,mag_column[1])     
    left_pipette.aspirate(300,reagent_column[9])  
    left_pipette.dispense(300,mag_column[1])
    left_pipette.drop_tip()


    #Second wash step of Col1
    left_pipette.pick_up_tip()
    left_pipette.aspirate(300,mag_column[0])
    left_pipette.dispense(300,reagent_column[5])     
    left_pipette.aspirate(300,mag_column[0])
    left_pipette.dispense(300,reagent_column[5])  #bin waste
    left_pipette.drop_tip()
    left_pipette.pick_up_tip()
    left_pipette.aspirate(300,reagent_column[10])  #adding wash buffer to col 1
    left_pipette.dispense(300,mag_column[0])     
    left_pipette.aspirate(300,reagent_column[10])  
    left_pipette.dispense(300,mag_column[0])
    left_pipette.drop_tip()


    #Second wash step of col2
    left_pipette.pick_up_tip()
    left_pipette.aspirate(300,mag_column[1])
    left_pipette.dispense(300,reagent_column[4])     
    left_pipette.aspirate(300,mag_column[0])
    left_pipette.dispense(300,reagent_column[4])  #bin waste
    left_pipette.drop_tip()
    left_pipette.pick_up_tip()
    left_pipette.aspirate(300,reagent_column[9])  #adding wash buffer to col 2
    left_pipette.dispense(300,mag_column[1])     
    left_pipette.aspirate(300,reagent_column[9])  
    left_pipette.dispense(300,mag_column[1])
    left_pipette.drop_tip()


    #Removing col1 final wash buffer
    left_pipette.well_bottom_clearance.dispense = 3
    left_pipette.well_bottom_clearance.aspirate = 0.2

    left_pipette.pick_up_tip()
    left_pipette.aspirate(300,mag_column[0])
    left_pipette.dispense(300,reagent_column[3])     
    left_pipette.aspirate(300,mag_column[0])
    left_pipette.dispense(300,reagent_column[3])  #bin waste
    left_pipette.drop_tip()

    #removing col2 final wash buffer
    protocol.delay(seconds = 30)
    left_pipette.pick_up_tip()
    left_pipette.aspirate(300,mag_column[1])
    left_pipette.dispense(300,reagent_column[2])     
    left_pipette.aspirate(300,mag_column[1])
    left_pipette.dispense(300,reagent_column[2])  #bin waste
    left_pipette.drop_tip()

    left_pipette.well_bottom_clearance.dispense = 1
    left_pipette.well_bottom_clearance.dispense = 1


    protocol.delay(minutes = 15)

    
    #add elution buffer to col1
    left_pipette.pick_up_tip()
    mag_deck.disengage()
    left_pipette.aspirate(40,reagent_column[0])
    left_pipette.dispense(40,mag_column[0])
    left_pipette.mix(5,30)
    left_pipette.drop_tip()

    left_pipette.pick_up_tip()
    left_pipette.aspirate(40,reagent_column[0])
    left_pipette.dispense(40,mag_column[1])
    left_pipette.mix(5,30)
    left_pipette.drop_tip()

    protocol.delay(minutes = 2.5)
    mag_deck.engage()
    protocol.delay(minutes = 2)    #Now DNA is purified, can aspirate and dispense anyway

