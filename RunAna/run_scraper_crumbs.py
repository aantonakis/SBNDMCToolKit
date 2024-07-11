import os
import sys

input_file = sys.argv[1]

output_file = sys.argv[2]

os.system("lar -c run_MCScraper.fcl -s "+input_file+" -o "+output_file)



