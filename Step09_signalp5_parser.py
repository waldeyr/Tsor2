'''
What does this program do?
	It parses signalp results and save them into a file temp02, which will be merged with other results

Input:
    Tsor2_signalp_summary.signalp5

Output:
	temp02.tsv

@author: Waldeyr Mendes Cordeiro da Silva, June 2020
'''


import pandas as pd
import re

signalp5Results = pd.read_csv('Tsor2_signalp_summary.signalp5', sep="\t")

head = ['contig_id', 'SignalP_Prediction', 'SignalP_SP', 'SignalP_Other', 'SignalP_Cleavage_Site_Position',
        'SignalP_Contig_Frame']
#print(signalp5Results)
# function to return signalp resuls as a formated string
def formatResult(value, frame):
	if (str(value['SignalP_Prediction']) != "OTHER"):
		positive_probability = str(value['SignalP_SP']).replace(",", ".")
		# negative_probability = str(value['SignalP_Other'])
		position_site = re.search(r"(\d+-\d+)", str(value['SignalP_Cleavage_Site_Position'])).group()
		cleavage_site = re.search(
			r"([A|R|N|D|C|E|Q|G|H|I|L|K|M|F|P|S|T|W|Y|V|\*]+-[A|R|N|D|C|E|Q|G|H|I|L|K|M|F|P|S|T|W|Y|V|\*]+)",
			str(value['SignalP_Cleavage_Site_Position'])).group()
		return f"Frame:{frame}; Positive Probability:{positive_probability}; Position: {position_site}; Cleavage site: {cleavage_site}. "
	else:
		#return f"Frame:{frame}; None. "
		return f""

SignalP_contig_id = []
SignalP_prediction = []
last_contig_id = str(signalp5Results['contig_id'][0])[:-2].strip()
resultFrames=resultFrame1=resultFrame2=resultFrame3=resultFrame4=resultFrame5=resultFrame6 = ""
# for each Signalp5 result
for index, value in signalp5Results.iterrows():
	current_contig_id = str(value['contig_id'])[:-2].strip()
	currentFrame = str(value['contig_id'])[-1].strip()
	if (currentFrame == "1"):
		resultFrame1 = formatResult(value, currentFrame)
	if (currentFrame == "2"):
		resultFrame2 = formatResult(value, currentFrame)
	if (currentFrame == "3"):
		resultFrame3 = formatResult(value, currentFrame)
	if (currentFrame == "4"):
		resultFrame4 = formatResult(value, currentFrame)
	if (currentFrame == "5"):
		resultFrame5 = formatResult(value, currentFrame)
	if (currentFrame == "6"):
		resultFrame6 = formatResult(value, currentFrame)
	if( (index+1)%6 == 0):
		resultFrames = resultFrame1+resultFrame2+resultFrame3+resultFrame4+resultFrame5+resultFrame6
		SignalP_contig_id.append(current_contig_id.strip())
		SignalP_prediction.append(resultFrames.strip())
		resultFrames = ""
	last_contig_id = current_contig_id

data = {
	"contig_id": SignalP_contig_id,
	"SignalP_prediction": SignalP_prediction
}

signalp5Adjusted = pd.DataFrame(data)
print(signalp5Adjusted.shape)
signalp5Adjusted.reset_index(drop=True, inplace=True)
signalp5Adjusted.set_index(['contig_id'], inplace = True)
signalp5Adjusted.to_csv("temp02.tsv", sep='\t', encoding='utf-8')
