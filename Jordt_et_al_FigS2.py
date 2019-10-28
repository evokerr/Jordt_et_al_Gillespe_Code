# This code uses the Gillespe algorithm to simulate governing
# differential equations in Jordt et al. to produce the data
# for Figure S2

import numpy as np
import csv

# Model parameters, state variables, running conditions and output files
mu=0.7
cost=0.2
K=0.004
Y=0.000002
lambda_v=0.05
tao=0.00001
agents=["A","Aa","Ab","Aab","B","Ba","Bb","Bab","R"]
agent_of_index={0:"A",1:"Aa",2:"Ab",3:"Aab",4:"B",5:"Ba",6:"Bb",7:"Bab", 8:"R"}
dilution_factor=0.02
batch_resources=0.02
steps=100000
period=8
stop_time=period*4
print_period=.1
init_cells=200
N_reps=20
csvfile = '/Users/Jordt_et_al_FigS2.csv'

# Parameters to vary
PS=[0.0,0.99] # corresponds to "p" values
PM=[0.0] # corresponds to "q" values


# Functions

# CalculateGrowths computes the weights for the growth of each
# bacterial type given the current population configuration
def CalculateGrowths(n, p):
	weight=[0 for i in range(8)]
	weight[0]=(
	((p["mu_A"]*n["R"])/(p["K"]+n["R"]))*n["A"]+
	((p["lambda"]*p["mu_Aa"]*n["R"])/(p["K"]+n["R"]))*n["Aa"]+
	((p["lambda"]*p["mu_Ab"]*n["R"])/(p["K"]+n["R"]))*n["Ab"]+
	((p["lambda"]*p["lambda"]*p["mu_Aab"]*n["R"])/(p["K"]+n["R"]))*n["Aab"]
	)
	weight[1]=(
	(((1.0-p["lambda"])*p["mu_Aa"]*n["R"])/(p["K"]+n["R"]))*n["Aa"]+
	(((1.0-p["lambda"])*p["lambda"]*p["mu_Aab"]*n["R"])/(p["K"]+n["R"]))*n["Aab"]
	)
	weight[2]=(
	(((1.0-p["lambda"])*p["mu_Ab"]*n["R"])/(p["K"]+n["R"]))*n["Ab"]+
	(((1.0-p["lambda"])*p["lambda"]*p["mu_Aab"]*n["R"])/(p["K"]+n["R"]))*n["Aab"]
	)
	weight[3]=(
	(((1.0-p["lambda"])*(1.0-p["lambda"])*p["mu_Aab"]*n["R"])/(p["K"]+n["R"]))*n["Aab"]
	)
	weight[4]=(
	((p["mu_B"]*n["R"])/(p["K"]+n["R"]))*n["B"]+
	((p["lambda"]*p["mu_Ba"]*n["R"])/(p["K"]+n["R"]))*n["Ba"]+
	((p["lambda"]*p["mu_Bb"]*n["R"])/(p["K"]+n["R"]))*n["Bb"]+
	((p["lambda"]*p["lambda"]*p["mu_Bab"]*n["R"])/(p["K"]+n["R"]))*n["Bab"]
	)
	weight[5]=(
	(((1.0-p["lambda"])*p["mu_Ba"]*n["R"])/(p["K"]+n["R"]))*n["Ba"]+
	(((1.0-p["lambda"])*p["lambda"]*p["mu_Bab"]*n["R"])/(p["K"]+n["R"]))*n["Bab"]
	)
	weight[6]=(
	(((1.0-p["lambda"])*p["mu_Bb"]*n["R"])/(p["K"]+n["R"]))*n["Bb"]+
	(((1.0-p["lambda"])*p["lambda"]*p["mu_Bab"]*n["R"])/(p["K"]+n["R"]))*n["Bab"]
	)
	weight[7]=(
	(((1.0-p["lambda"])*(1.0-p["lambda"])*p["mu_Bab"]*n["R"])/(p["K"]+n["R"]))*n["Bab"]
	)
	return(weight)

# CalculateTransfers computes the weights for the transfer of plasmids
# to every bacterial type given the current population configuration
def CalculateTransfers(n, p):
	weight=[0 for i in range(8)]
	weight[0]=(p["tao"]/2)*(2*n["Aa"]+n["Aab"]+2*n["Ba"]+n["Bab"])*n["A"]
	weight[1]=(p["tao"]/2)*(2*n["Ab"]+n["Aab"]+2*n["Bb"]+n["Bab"])*n["A"]
	weight[2]=(p["tao"]/2)*(2*n["Ab"]+n["Aab"]+2*n["Bb"]+n["Bab"])*n["Aa"]
	weight[3]=(p["tao"]/2)*(2*n["Aa"]+n["Aab"]+2*n["Ba"]+n["Bab"])*n["Ab"]
	weight[4]=(p["tao"]/2)*(2*n["Aa"]+n["Aab"]+2*n["Ba"]+n["Bab"])*n["B"]
	weight[5]=(p["tao"]/2)*(2*n["Ab"]+n["Aab"]+2*n["Bb"]+n["Bab"])*n["B"]
	weight[6]=(p["tao"]/2)*(2*n["Ab"]+n["Aab"]+2*n["Bb"]+n["Bab"])*n["Ba"]
	weight[7]=(p["tao"]/2)*(2*n["Aa"]+n["Aab"]+2*n["Ba"]+n["Bab"])*n["Bb"]
	return(weight)


# Main Loop

with open(csvfile, "w") as output2:
	writer = csv.writer(output2, lineterminator='\n')
	writer.writerows([["p"]+["q"]+["rep"]+["time"]+agents])
	x=-1
	for psv in PS:
		x+=1
		ps=psv
		y=-1
		for pmv in PM:
			y+=1
			pm=pmv
			value_of_parameter={"mu_A": mu, "mu_Aa": mu-(1.0-ps)*cost, "mu_Ab": mu-(1.0-pm)*cost, 
			"mu_Aab": mu-(1.0-ps)*cost-(1.0-pm)*cost, "mu_B": mu, "mu_Ba": mu-(1.0-pm)*cost, 
			"mu_Bb": mu-(1.0-ps)*cost, "mu_Bab": mu-(1.0-ps)*cost-(1.0-pm)*cost, 
			"K": K, "Y": Y, "lambda": lambda_v, "tao": tao}
			for r in range(N_reps):
				number_of_agent={"A":0 ,"Aa":init_cells/2 ,"Ab":0 ,"Aab":0 ,"B":0 ,"Ba":0 ,"Bb":init_cells/2 ,"Bab":0 ,"R":batch_resources}
				t=0
				p=1
				pp=1
				while t < stop_time:
					growth=CalculateGrowths(number_of_agent, value_of_parameter)
					transfer=CalculateTransfers(number_of_agent, value_of_parameter)
					sg=sum(growth)
					st=sum(transfer)
					t=t+np.random.exponential(1/(sg+st))
					if np.random.binomial(1,sg/(sg+st))==1:
						number_of_agent["R"]=number_of_agent["R"]-(value_of_parameter["Y"])
						who_grows=np.random.multinomial(1,[growth[j]/sg for j in range(8)]).tolist()
						for j in range(8):
							if who_grows[j]==1:
								index=j
						number_of_agent[agent_of_index[index]]+=1			
					else:
						who_tranfers=np.random.multinomial(1,[transfer[j]/st for j in range(8)]).tolist()
						for j in range(8):
							if who_tranfers[j]==1:
								index=j
						if index==0:
							number_of_agent["A"]=number_of_agent["A"]-1
							number_of_agent["Aa"]=number_of_agent["Aa"]+1
						if index==1:
							number_of_agent["A"]=number_of_agent["A"]-1
							number_of_agent["Ab"]=number_of_agent["Ab"]+1		
						if index==2:
							number_of_agent["Aa"]=number_of_agent["Aa"]-1
							number_of_agent["Aab"]=number_of_agent["Aab"]+1			
						if index==3:
							number_of_agent["Ab"]=number_of_agent["Ab"]-1
							number_of_agent["Aab"]=number_of_agent["Aab"]+1	
						if index==4:
							number_of_agent["B"]=number_of_agent["B"]-1
							number_of_agent["Ba"]=number_of_agent["Ba"]+1
						if index==5:
							number_of_agent["B"]=number_of_agent["B"]-1
							number_of_agent["Bb"]=number_of_agent["Bb"]+1
						if index==6:
							number_of_agent["Ba"]=number_of_agent["Ba"]-1
							number_of_agent["Bab"]=number_of_agent["Bab"]+1		
						if index==7:
							number_of_agent["Bb"]=number_of_agent["Bb"]-1
							number_of_agent["Bab"]=number_of_agent["Bab"]+1
					if t > p*period:
						for j in range(8):
							number_of_agent[agent_of_index[j]]=np.random.binomial(number_of_agent[agent_of_index[j]],dilution_factor)
						number_of_agent["R"]=batch_resources
						p+=1
					if t > pp*print_period:
						writer.writerows([[ps]+[pm]+[r]+[t]+[number_of_agent[agent_of_index[j]] for j in range(9)]])
						pp+=1
				print(r)
			print([ps,pm])
