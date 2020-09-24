/*********************************************
 * OPL 12.6.2.0 Data
 * Author: XYZ
 * Creation Date: Jan 23, 2018 at 10:09:07 AM
 *********************************************/

/**** Raw Data Tuple ****/
tuple Raw_Data{
	string rownames;
	int FY;
	int QoY;
	int MoY;
	string dmu;// thus is the dmu
	string CallCenterType;
	string ServiceAreaDesc;
	string CallCenterDesc;
	string FunctionGroup;
	string ParentFunctionDesc;
	float TotalOfferedCalls;
	float TotalCallsAnswered;
	float TotalHandleTime;
	float TotalTransferredCalls;
	float TotalTachTransfers;
	float Repeats2Hour;
	float Repeats3Day;
	float NTF;
	float Attrition;
	float Reliably;
	float Occupancy;
	float TotalTalkTime;
	float TotalHoldTime;
	float TotalAfterCallWorkTime;
};

{Raw_Data} RawData = ...; 

/***** Scenario Parameters ****/
string scenario = "Test_Nov17_2";
string fiscYear;
string month;
string quarter;
string Pfunction;
/******************************/

execute{
	var temp = new Date();
	//scenario = temp ;
	for(var i in RawData){
		fiscYear = i.FY;
		month = i.MoY;
		quarter = i.QoY;
		Pfunction = i.ParentFunctionDesc; 
	}
}

ordered {string} DMUs = {j.dmu | j in RawData};

int currIter =...;
int TotalDMUs;

execute{
	TotalDMUs = DMUs.size;
}
string currDMUs[1..TotalDMUs];
string currDMU = currDMUs[currIter];

// In case runningv this code for a single Call Center rather than all of them
// below line needs to get acticated to manually overwrite the DMU parameter
//string currDMU;//=...;// ="XYZ";

execute POP_currDMUs {
	var count=1;
	for(var x in DMUs){
		currDMUs[count]=x;
		count++;	
	}
}

// currDMU (Decision Making Unit- to be consistent w/ the literature) is the Call Center
// being rated, we iterate through all Call Centers and solve this model for each CenterID
// Below paramters needs to get adjusted by scenario as the scale matters
float epsilon_tech = 0.00000002;  //a small number to keep the Obj positive
float epsilon_care = 0.000000001;
float epsilon_care_quart = 0.0000000005;
float epsilon_tech_jan2018 = 0.0000000101; 

// DEA Efficiency and Linear Combination Variables 
dvar float theta; //efficiency rating it will be between 0 and 1
dvar float+ lambda[DMUs];
// this is the linear combination of other Call Centerss, or virtual Call Center
//, used to compare to the current Call Center, or currDMU in model

// DEA inputs are below variables
dvar float+ TotalHandleTimeSlack; // Labor Hour
dvar float+ TotalOfferedCallsSlack; // Total Traffic Offered i.e. Calls Offered

// DEA outputs are below variables
dvar float+ TotalCallsAnsweredSurplus;
dvar float+ TotalTransferredCallsSurplus;
dvar float+ TotalTachTransfersSurplus;
dvar float+ Repeats2HourSurplus;
dvar float+ Repeats3DaySurplus;
dvar float+ NTFSurplus;
dvar float+ AttritionSurplus;
dvar float+ ReliablySurplus;
dvar float+ OccupancySurplus;
dvar float+ TotalTalkTimeSurplus;
dvar float+ TotalHoldTimeSurplus;
dvar float+ TotalAfterCallWorkTimeSurplus;

dexpr float TotalSurplus = TotalHandleTimeSlack + TotalOfferedCallsSlack 
							+ TotalCallsAnsweredSurplus + TotalTransferredCallsSurplus 
							+ TotalTachTransfersSurplus + Repeats2HourSurplus + Repeats3DaySurplus
							+ NTFSurplus + AttritionSurplus + ReliablySurplus + OccupancySurplus
							+ TotalTalkTimeSurplus + TotalHoldTimeSurplus + TotalAfterCallWorkTimeSurplus;

/***** This is the Dual Formulation of the DEA Model ****/
minimize theta - epsilon_care*(TotalSurplus);

subject to{
	// we should have inequality constraints but w/ dummy slack and surplus DVs we can easily track the slack/surplus in each constraint
	ct_1:
		sum(d in RawData: d.dmu==currDMU)(d.TotalHandleTime)*theta == sum(d in RawData) d.TotalHandleTime*lambda[d.dmu] + TotalHandleTimeSlack;
	
	/*ct_2:	
		sum(d in RawData: d.dmu==currDMU)(d.TotalOfferedCalls)*theta == sum(d in RawData) d.TotalOfferedCalls*lambda[d.dmu] + TotalOfferedCallsSlack;
		*/	 
	ct_3:
		sum(d in RawData: d.dmu==currDMU) d.TotalCallsAnswered == sum(d in RawData) d.TotalCallsAnswered*lambda[d.dmu] - TotalCallsAnsweredSurplus;
	
	ct_4:
		sum(d in RawData: d.dmu==currDMU) d.TotalTransferredCalls == sum(d in RawData) d.TotalTransferredCalls*lambda[d.dmu] - TotalTransferredCallsSurplus;
	
	ct_5:
		sum(d in RawData: d.dmu==currDMU) d.TotalTachTransfers == sum(d in RawData) d.TotalTachTransfers*lambda[d.dmu] - TotalTachTransfersSurplus;
	
	ct_6:
		sum(d in RawData: d.dmu==currDMU) d.Repeats2Hour  == sum(d in RawData) d.Repeats2Hour*lambda[d.dmu] - Repeats2HourSurplus;
	
	ct_7:
		sum(d in RawData: d.dmu==currDMU) d.Repeats3Day  == sum(d in RawData) d.Repeats3Day*lambda[d.dmu] - Repeats3DaySurplus;
	
	ct_8:
		sum(d in RawData: d.dmu==currDMU) d.NTF  == sum(d in RawData) d.NTF*lambda[d.dmu] - NTFSurplus;		
	
	ct_9:		
		sum(d in RawData: d.dmu==currDMU) d.Attrition  == sum(d in RawData) d.Attrition*lambda[d.dmu] - AttritionSurplus;
			
	ct_10:
		sum(d in RawData: d.dmu==currDMU) d.Reliably == sum(d in RawData) d.Reliably*lambda[d.dmu] - ReliablySurplus;
			
	ct_11:
		sum(d in RawData: d.dmu==currDMU) d.Occupancy  == sum(d in RawData) d.Occupancy*lambda[d.dmu] - OccupancySurplus;
		
	/*ct_12:
		sum(d in RawData: d.dmu==currDMU)(d.TotalTalkTime) == sum(d in RawData) d.TotalTalkTime*lambda[d.dmu] - TotalTalkTimeSurplus;
		
	ct_13:
		sum(d in RawData: d.dmu==currDMU)(d.TotalHoldTime) == sum(d in RawData) d.TotalHoldTime *lambda[d.dmu] - TotalHoldTimeSurplus;	*/
			 	
	ct_14:
		sum(d in RawData: d.dmu==currDMU)(d.TotalAfterCallWorkTime) == sum(d in RawData) d.TotalAfterCallWorkTime*lambda[d.dmu] - TotalAfterCallWorkTimeSurplus;	 
		
	ct_15: // this constraint forces model to only create virtual Call Centers from linear combinations of Call Centers which are of similar size
		sum(d in DMUs)lambda[d]== 1; //this is appended to turn CCR model into BBC model
};

tuple ThetaSolution {
	string dmu; //current decision making unit, or call center
	string ParentFunctionDesc;
	int fiscal_year;
	int month;
	string scenario;
	float theta; //efficiency value
}
{ThetaSolution} thetaSol;

tuple LambdaSolution {
	string dmu; //current decision making unit, or call center
	string ParentFunctionDesc;
	string scenario;
	int fiscal_year;
	int month;
	string callcenter; //call center used to create virtual call center
	float lambda; //proportion of call center used in virtual call center
}
{LambdaSolution} lambdaSol;

execute FormatSolution {
	thetaSol.add(currDMU, Pfunction, fiscYear, month, scenario, theta)	
	for(var d in DMUs){
		if(lambda[d]>0){
			lambdaSol.add(currDMU, Pfunction, scenario, fiscYear, month, d, lambda[d])	
		}
	}
}

tuple effiSolutions {
	string dmu;
	string ParentFunctionDesc;
	string scenario;
	int fiscal_year;
	int month;
	float TotalHandleTimeSlack;
	float TotalOfferedCallsSlack;
	float TotalCallsAnsweredSurplus;
	float TotalTransferredCallsSurplus;
	float TotalTachTransfersSurplus;
	float Repeats2HourSurplus;
	float Repeats3DaySurplus;
	float NTFSurplus;
	float AttritionSurplus;
	float ReliablySurplus;
	float OccupancySurplus;
	float TotalTalkTimeSurplus;
	float TotalHoldTimeSurplus;
	float TotalAfterCallWorkTimeSurplus;
}
{effiSolutions} EfficiencySolution;

execute{
	EfficiencySolution.add(currDMU, Pfunction, scenario, fiscYear, month
					, TotalHandleTimeSlack.solutionValue, TotalOfferedCallsSlack.solutionValue
					, TotalCallsAnsweredSurplus.solutionValue
					, TotalTransferredCallsSurplus.solutionValue
					, TotalTachTransfersSurplus.solutionValue, Repeats2HourSurplus.solutionValue
					, Repeats3DaySurplus.solutionValue, NTFSurplus.solutionValue
					, AttritionSurplus.solutionValue, ReliablySurplus.solutionValue
					, OccupancySurplus.solutionValue, TotalTalkTimeSurplus.solutionValue
					, TotalHoldTimeSurplus.solutionValue, TotalAfterCallWorkTimeSurplus.solutionValue);
}

main {
	thisOplModel.generate();
	var constantData =  new IloOplDataElements();
	// making a variable to hold the tuple set "RawData" which holds the data input table
	// "RawData does not change, i.e., is constant, and I dont want to read it in at each iteration"
	constantData.RawData = thisOplModel.RawData;
	
	// similar idea for allDMUs as for "RawData"
	var allDMUs = thisOplModel.DMUs;
	
	// Path to the source .mod file
	var source = new IloOplModelSource("C:\\Users\\abdolho\\Desktop\\CallCenter\\DEA\\SPC_DEA\\SPC_DEA.mod");
	
	//This file path to .dat file I use to store solution on SQL DB 
	var solData = new IloOplDataSource("C:\\Users\\abdolho\\Desktop\\CallCenter\\DEA\\SPC_DEA\\SPC_DEA_Write2DB.dat");
	
	//Create the model defination
	var def = new IloOplModelDefinition(source);
	
	//Variable to hold all the lambda solutions and initialize to all zeros
	var allLambdas = new Array();
	
	for(var x in allDMUs){
		allLambdas[x]= new Array();
		for(var y in allDMUs){
			allLambdas[x][y]=0;	
		}
	}
	// Dont need this part necessarily 
	var effFrontiers = new Array();
	for(var x in allDMUs){
		effFrontiers[x] = new Array();
		for(var y = 1 ; y <= 10 ; y++){
			effFrontiers[x][y]=0;	
		}
	}
	// This is the main loop to solve the LP iteratively 
	for(var i = 1 ; i <= thisOplModel.TotalDMUs ; i++){
		// have to solve the model once for each Call Center	
		var cplex = new IloCplex();
		var opl = new IloOplModel(def,cplex);
		var myData = new IloOplDataElements(); 
		myData.currIter = i;
		//myData.currDMU=thisOplModel.currDMUs[i];
		opl.addDataSource(constantData);
		opl.addDataSource(myData);
		opl.addDataSource(solData); //had to add this so that I can call DBUpdate() from main block
		opl.generate();
		cplex.solve();
		writeln(opl.currDMUs[i]," > ",opl.theta.solutionValue);	
		for(var d in allDMUs){
			allLambdas[opl.currDMUs[i]][d] = opl.lambda[d].solutionValue;
		}
		
		writeln("*************************************");
		writeln("Total Handle Time Slack -> ", opl.TotalHandleTimeSlack.solutionValue);
		writeln("Total Offered Calls Slack -> ", opl.TotalOfferedCallsSlack.solutionValue);
		writeln("Total Calls Answered Surplus -> ", opl.TotalCallsAnsweredSurplus.solutionValue);
		writeln("Total Transferred Calls Surplus -> ", opl.TotalTransferredCallsSurplus.solutionValue);
		writeln("Total Tech Transfers Surplus -> ", opl.TotalTachTransfersSurplus.solutionValue);
		writeln("Repeats 2Hour Surplus -> ", opl.Repeats2HourSurplus.solutionValue);
		writeln("Repeats 3Day Surplus -> ", opl.Repeats3DaySurplus.solutionValue);
		writeln("NTF Surplus -> ", opl.NTFSurplus.solutionValue);
		writeln("Attrition Surplus -> ", opl.AttritionSurplus.solutionValue);
		writeln("Reliably Surplus -> ", opl.ReliablySurplus.solutionValue);
		writeln("Occupancy Surplus -> ", opl.OccupancySurplus.solutionValue);
		writeln("Total Talk Time Surplus -> ", opl.TotalTalkTimeSurplus.solutionValue);
		writeln("Total Hold Time Surplus -> ", opl.TotalHoldTimeSurplus.solutionValue);
		writeln("Total After Call Work Time Surplus -> ", opl.TotalAfterCallWorkTimeSurplus.solutionValue);
		writeln("*************************************");
		opl.postProcess();
		opl.end();
		cplex.end(); 
		myData.end();
	}
	writeln("printing lambdas...");
	for(var d in allDMUs){
		writeln("Lambdas_for_",d);
		for(var e in allDMUs){
			if(allLambdas[d][e] > 0){	
				write("lambda_of>", e,">");
				writeln(allLambdas[d][e]);	
 			}				
		}
	}
	source.end();
	def.end()
	constantData.end(); 
} 