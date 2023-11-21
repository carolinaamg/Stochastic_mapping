import sys
import os
import re

input_folder = sys.argv[1] ##Folder where the output files of 1_stochastic_mapping_random.R are
output_file = open(sys.argv[2],"w") ##Output file where distribution will go

output_file.write("Analysis" + "\t" + "Node" + "\t" + "Replicate" + "\t" + "Clusters_gained_or_lost" + "\n")
for file in os.listdir(input_folder):
    if "gains_stochastic_SM.tsv" in file:
        gains_stochastic = os.path.join(input_folder, file)
        losses_stochastic = gains_stochastic.replace("gains","losses")
        gains_simulation = gains_stochastic.replace("stochastic","simulation")
        losses_simulation = gains_stochastic.replace("gains","losses")

        ##Doing gains stochastic
        with open(gains_stochastic) as gains_stochastic:
            for line in gains_stochastic:
                line = line.rstrip()

                if line.startswith("Node_"):
                    tabs = line.split("\t")
                    node_name = tabs [0]
                    replicates = tabs [1:len(tabs)]

                    count = 0
                    for x in replicates:
                        count = count + 1
                        if x == "character(0)":
                            result = "0"
                            output_file.write("Gains_stochastic" + "\t" + node_name + "\t" + str(count) + "\t" + result + "\n")
                        else:
                            x = x.replace("c(","")
                            x = x.replace('"','')
                            x = x.replace(")","")
                            x = x.split(", ")
                            result = ";".join(x)
                            output_file.write("Gains_stochastic" + "\t" + node_name + "\t" + str(count) + "\t" + result + "\n")
                else:
                    pass
        gains_stochastic.close()

        ##Doing losses stochastic
        with open(losses_stochastic) as losses_stochastic:
            for line in losses_stochastic:
                line = line.rstrip()

                if line.startswith("Node_"):
                    tabs = line.split("\t")
                    node_name = tabs [0]
                    replicates = tabs [1:len(tabs)]

                    count = 0
                    for x in replicates:
                        count = count + 1
                        if x == "character(0)":
                            result = "0"
                            output_file.write("Losses_stochastic" + "\t" + node_name + "\t" + str(count) + "\t" + result + "\n")
                        else:
                            x = x.replace("c(","")
                            x = x.replace('"','')
                            x = x.replace(")","")
                            x = x.split(", ")
                            result = ";".join(x)
                            output_file.write("Losses_stochastic" + "\t" + node_name + "\t" + str(count) + "\t" + result + "\n")
                else:
                    pass
        losses_stochastic.close()

        ##Doing gains simulations
        with open(gains_simulation) as gains_simulation:
            for line in gains_simulation:
                line = line.rstrip()

                if line.startswith("Node_"):
                    tabs = line.split("\t")
                    node_name = tabs [0]
                    replicates = tabs [1:len(tabs)]

                    count = 0
                    for x in replicates:
                        count = count + 1
                        if x == "character(0)":
                            result = "0"
                            output_file.write("Gains_simulation" + "\t" + node_name + "\t" + str(count) + "\t" + result + "\n")
                        else:
                            x = x.replace("c(","")
                            x = x.replace('"','')
                            x = x.replace(")","")
                            x = x.split(", ")
                            result = ";".join(x)
                            output_file.write( "Gains_simulation" + "\t" + node_name + "\t" + str(count) + "\t" + result + "\n")
                else:
                    pass
        gains_simulation.close()

        ##Doing losses simulations
        with open(losses_simulation) as losses_simulation:
            for line in losses_simulation:
                line = line.rstrip()

                if line.startswith("Node_"):
                    tabs = line.split("\t")
                    node_name = tabs [0]
                    replicates = tabs [1:len(tabs)]

                    count = 0
                    for x in replicates:
                        count = count + 1
                        if x == "character(0)":
                            result = "0"
                            output_file.write("Losses_simulation" + "\t" + node_name + "\t" + str(count) + "\t" + result + "\n") 
                        else:
                            x = x.replace("c(","")
                            x = x.replace('"','')
                            x = x.replace(")","")
                            x = x.split(", ")
                            result = ";".join(x)
                            output_file.write("Losses_simulation" + "\t" + node_name + "\t" + str(count) + "\t" + result + "\n")
                else:
                    pass
        losses_stochastic.close()
output_file.close()