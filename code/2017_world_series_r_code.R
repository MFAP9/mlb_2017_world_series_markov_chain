library(markovchain)

setwd("~/Northwestern University/Summer 2018 Quarter/MSDS 456/Assignments/Assignment 3/play_files")
# READ IN THE CSV FILES
fn<-list.files("~/Northwestern University/Summer 2018 Quarter/MSDS 456/Assignments/Assignment 3/play_files")

for (i in 1:30){
  fname<-fn[i]
  fname<-substr(fname,start = 5, stop = 7)
  assign(fname,read.csv(fn[i],header = FALSE))
  }
fullset<-rbind(ANA,ARI,ATL,BAL,BOS,CHA,CHN,CIN,CLE,COL,DET,HOU,KCA,LAN,MIA,MIL,MIN,NYA,NYN,OAK,PHI,PIT,SDN,SEA,SFN,SLN,TBA,TEX,TOR,WAS)

names(fullset)<-c("Game_ID","Visiting_Team","Inning","Batting_Team","Outs","Balls","Strikes","Vistor_Score","Home_Score","Batter","Batter_Hand","Pitcher","Pitcher_Hand","First_Runner","Second_Runner","Third_Runner","Event_Text","Leadoff_Flag","PinchHit_Flag","Defensive_Position","Lineup_Position","Event_Type","Batter_Event_Flag","AB_Flag","Hit_Value","SH_Flag","SF_Flag","Outs_On_Play","RBIs_On_Play","Wild_Pitch_Flag","Passed_Ball_Flag","Num_Errors","Batter_Dest","Runner_On_1st_Dest","Runner_On_2nd_Dest","Runner_On_3rd_Dest")

#ADD INDICATORS FOR RUNNERS AT EACH BASE
fullset["First"]<-0
fullset["Second"]<-0
fullset["Third"]<-0
fullset["Runs_Scored"]<-0
fullset["Result_Outs"]<-fullset$Outs+fullset$Outs_On_Play

for (i in 1:nrow(fullset)){
  if (fullset$First_Runner[i]!=""){
    fullset$First[i]<-1
  }
  if (fullset$Second_Runner[i]!=""){
    fullset$Second[i]<-1
  }
  if (fullset$Third_Runner[i]!=""){
    fullset$Third[i]<-1
  }
  if (fullset$Batter_Dest[i]>=4){
    fullset$Runs_Scored[i]<-fullset$Runs_Scored[i]+1
  }
  if (fullset$Runner_On_1st_Dest[i]>=4){
    fullset$Runs_Scored[i]<-fullset$Runs_Scored[i]+1
  }
  if (fullset$Runner_On_2nd_Dest[i]>=4){
    fullset$Runs_Scored[i]<-fullset$Runs_Scored[i]+1
  }
  if (fullset$Runner_On_3rd_Dest[i]>=4){
    fullset$Runs_Scored[i]<-fullset$Runs_Scored[i]+1
  }
}

#CREATE A VARIABLE TO INDICATE THE STARTING STATE OF EACH PLAY AND THE NEXT STATE
#STATES ARE OUTS,FIRST,SECOND,THIRD SO 1010 WOULD BE ONE OUT RUNNER ON SECOND ; 2101 WOULD BE TWO OUTS WITH RUNNERS ON FIRST AND THIRD; 3000 INDICATES THREE OUTS HAVE BEEN REACHED
fullset["State"]<-0
library(cat)
fullset$State<-paste(fullset$Outs,fullset$First,fullset$Second,fullset$Third,sep = "")

fullset["Next_State"]<-0
for (i in 1:nrow(fullset)){
  j<-i+1
  if (fullset$Result_Outs[i]==3){
    fullset$Next_State[i]<-3000
  }
  else{
    fullset$Next_State[i]<-fullset$State[j]
  }
}

#MAKES TWO WAY FREQUENCY TABLE OF TRANSITION STATES
library(gmodels)
twotab<-as.data.frame(table(fullset$State,fullset$Next_State))
names(twotab)<-c("State","Next_State","Freq")
aggtable<-aggregate(twotab$Freq,by=list(twotab$State),FUN=sum)
names(aggtable)<-c("State","Total")
twotab<-merge(twotab,aggtable,by="State")
twotab["Trans_Prob"]<-twotab$Freq/twotab$Total


sub<-subset(fullset,select = c("State","Next_State"))
twotab<-table(sub)
twotab<-table(sub)
twentyfive<-as.data.frame(c(100,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0))
names(twentyfive)<-c("3000")
twentyfive<-t(twentyfive)
twotab<-rbind(twotab,twentyfive)
twotab<-as.data.frame(twotab)

#CONVERTS TO PROBABILITY TRANSITION MATRIX AND THEN INTO A MARKOV CHAIN
twotab["Sum"]<-0
for (i in 1:25){
  twotab$Sum<-twotab$Sum+twotab[,i]
}
mlb.twotab<-twotab
for (i in 1:25){
  twotab[,i]<-twotab[,i]/twotab$Sum
}

twotab<-twotab[,1:25]
test<-subset(fullset,fullset$State=="2111"&fullset$Next_State=="0000")



states.list<-colnames(twotab)

twotab.m<-as.matrix(twotab)
MC<-new("markovchain",states=states.list,byrow=T,transitionMatrix=twotab.m,name="Baseball States")

#CALCULATES THE EXPECTED RUNS SCORED AS YOU MOVE FROM ONE STATE TO ANOTHER
exp.runs<-aggregate(fullset$Runs_Scored,by=list(fullset$State,fullset$Next_State),FUN=mean)
names(exp.runs)<-c("State","Next_State","Exp_Runs")

#LOOP THROUGH ALL POSSIBLE STATES TO CALCULATE HOW MANY RUNS YOU'D EXPECT TO SCORE BY THE END OF THE INNING GIVEN AN INITIAL STATE
all.states<-c("0000","0001","0010","0011","0100","0101","0110","0111","1000","1001","1010","1011","1100","1101","1110","1111","2000","2001","2010","2011","2100","2101","2110","2111","3000")
exp.table<-NA
for (x in 1:25){
  i<-all.states[x]
  tot.sim.runs<-0
  for (c in 1:10000){
    i<-all.states[x]
    vec.states<-c(i)
    iter<-0
    current.exp.runs<-0
    while(i!="3000"){

      mci<-as.data.frame(MC[i])
      names(mci)<-c("prob")
      mci["min"]<-0
      mci["max"]<-0
      mci["choice"]<-0
      mci$max[1]<-mci$prob[1]
      for (j in 2:nrow(mci)){
        k<-j-1
        mci$min[j]<-mci$max[k]
        mci$max[j]<-mci$min[j]+mci$prob[j]
      }
      mci$min<-mci$min*1000
      mci$max<-mci$max*1000
      #chance<-sample(1:1000,1)
      chance<-as.numeric(runif(1,min=0,max = 1000))
      for (a in 1:nrow(mci)){
        if ((chance>mci$min[a]) & (chance<=mci$max[a])){
          mci$choice[a]<-1
          n_state<-rownames(mci[a,])
        }
      }
      st.exp.runs<-subset(exp.runs,exp.runs$State==i & exp.runs$Next_State==n_state)
      current.exp.runs<-current.exp.runs+as.numeric(st.exp.runs$Exp_Runs)
      vec.states<-c(vec.states,n_state)
      i<-n_state
      iter<-iter+1
    }
    tot.sim.runs<-tot.sim.runs+current.exp.runs
  }
  avg.exp.runs<-tot.sim.runs/10000
  #print(all.states[x])
  #print(avg.exp.runs)
  row.add<-c(all.states[x],avg.exp.runs)
  exp.table<-rbind(exp.table,row.add)
}

exp.table.clean<-exp.table[-1,]
names(exp.table.clean)<-c("State","Exp_Runs")

######################################
######################################
######################################
#SUBSET ONLY AT-BATS BY THE ASTROS AND DODGERS
fullset["Home_Team"]<-substr(fullset$Game_ID,start = 1, stop = 3)
Dodgers<-subset(fullset,(fullset$Home_Team=="LAN" & fullset$Batting_Team==1) | (fullset$Visiting_Team=="LAN" & fullset$Batting_Team==0))
Astros<-subset(fullset,(fullset$Home_Team=="HOU" & fullset$Batting_Team==1) | (fullset$Visiting_Team=="HOU" & fullset$Batting_Team==0))

######################################
######################################
######################################
#DODGERS: THIS DOES THE SAME THING AS ABOVE BUT FOR ONLY DODGER AT-BATS
dod.twotab<-as.data.frame(table(Dodgers$State,Dodgers$Next_State))
names(dod.twotab)<-c("State","Next_State","Freq")
aggtable<-aggregate(dod.twotab$Freq,by=list(dod.twotab$State),FUN=sum)
names(aggtable)<-c("State","Total")
dod.twotab<-merge(dod.twotab,aggtable,by="State")
dod.twotab["Trans_Prob"]<-dod.twotab$Freq/dod.twotab$Total


sub<-subset(Dodgers,select = c("State","Next_State"))
dod.twotab<-table(sub)
dod.twotab<-table(sub)
twentyfive<-as.data.frame(c(100,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0))
names(twentyfive)<-c("3000")
twentyfive<-t(twentyfive)
dod.twotab<-rbind(dod.twotab,twentyfive)
dod.twotab<-as.data.frame(dod.twotab)

dod.twotab["Sum"]<-0
for (i in 1:25){
  dod.twotab$Sum<-dod.twotab$Sum+dod.twotab[,i]
}
for (i in 1:25){
  dod.twotab[,i]<-dod.twotab[,i]/dod.twotab$Sum
}

dod.twotab<-dod.twotab[,1:25]
test<-subset(Dodgers,Dodgers$State=="2111"&Dodgers$Next_State=="0000")



states.list<-colnames(dod.twotab)

dod.twotab.m<-as.matrix(dod.twotab)
dod.MC<-new("markovchain",states=states.list,byrow=T,transitionMatrix=dod.twotab.m,name="Baseball States - Dodgers")




dod.exp.runs<-aggregate(Dodgers$Runs_Scored,by=list(Dodgers$State,Dodgers$Next_State),FUN=mean)
names(dod.exp.runs)<-c("State","Next_State","Exp_Runs")

tempmerge<-merge(dod.exp.runs,exp.runs,by=c("State","Next_State"),all=T)
for (m in 1:nrow(tempmerge)){
  if (is.na(tempmerge$Exp_Runs.x[m])){
    tempmerge$Exp_Runs.x[m]<-tempmerge$Exp_Runs.y[m]
  }
}
tempmerge<-tempmerge[,1:3]
dod.exp.runs<-tempmerge
names(dod.exp.runs)<-c("State","Next_State","Exp_Runs")

#Loop all possible states
all.states<-c("0000","0001","0010","0011","0100","0101","0110","0111","1000","1001","1010","1011","1100","1101","1110","1111","2000","2001","2010","2011","2100","2101","2110","2111","3000")
dod.exp.table<-NA
for (x in 1:25){
  i<-all.states[x]
  tot.sim.runs<-0
  for (c in 1:10000){
    i<-all.states[x]
    vec.states<-c(i)
    iter<-0
    current.exp.runs<-0
    while(i!="3000"){
      
      mci<-as.data.frame(dod.MC[i])
      names(mci)<-c("prob")
      mci["min"]<-0
      mci["max"]<-0
      mci["choice"]<-0
      mci$max[1]<-mci$prob[1]
      for (j in 2:nrow(mci)){
        k<-j-1
        mci$min[j]<-mci$max[k]
        mci$max[j]<-mci$min[j]+mci$prob[j]
      }
      mci$min<-mci$min*1000
      mci$max<-mci$max*1000
      #chance<-sample(1:1000,1)
      chance<-as.numeric(runif(1,min=0,max = 1000))
      for (a in 1:nrow(mci)){
        if ((chance>mci$min[a]) & (chance<=mci$max[a])){
          mci$choice[a]<-1
          n_state<-rownames(mci[a,])
        }
      }
      st.exp.runs<-subset(dod.exp.runs,dod.exp.runs$State==i & dod.exp.runs$Next_State==n_state)
      current.exp.runs<-current.exp.runs+st.exp.runs$Exp_Runs
      vec.states<-c(vec.states,n_state)
      i<-n_state
      iter<-iter+1
    }
    tot.sim.runs<-tot.sim.runs+current.exp.runs
  }
  avg.exp.runs<-tot.sim.runs/10000
  #print(all.states[x])
  #print(tot.sim.runs)
  #print(avg.exp.runs)
  row.add<-c(all.states[x],avg.exp.runs)
  dod.exp.table<-rbind(dod.exp.table,row.add)
}

dod.exp.table.clean<-dod.exp.table[-1,]
names(dod.exp.table.clean)<-c("State","Exp_Runs")

######################################
######################################
######################################
#ASTROS: THIS DOES THE SAME THING AS ABOVE BUT JUST FOR ASTROS AT-BATS

ast.twotab<-as.data.frame(table(Astros$State,Astros$Next_State))
names(ast.twotab)<-c("State","Next_State","Freq")
aggtable<-aggregate(ast.twotab$Freq,by=list(ast.twotab$State),FUN=sum)
names(aggtable)<-c("State","Total")
ast.twotab<-merge(ast.twotab,aggtable,by="State")
ast.twotab["Trans_Prob"]<-ast.twotab$Freq/ast.twotab$Total


sub<-subset(Astros,select = c("State","Next_State"))
ast.twotab<-table(sub)
ast.twotab<-table(sub)
twentyfive<-as.data.frame(c(100,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0))
names(twentyfive)<-c("3000")
twentyfive<-t(twentyfive)
ast.twotab<-rbind(ast.twotab,twentyfive)
ast.twotab<-as.data.frame(ast.twotab)

ast.twotab["Sum"]<-0
for (i in 1:25){
  ast.twotab$Sum<-ast.twotab$Sum+ast.twotab[,i]
}
for (i in 1:25){
  ast.twotab[,i]<-ast.twotab[,i]/ast.twotab$Sum
}

ast.twotab<-ast.twotab[,1:25]
test<-subset(Astros,Astros$State=="2111"&Astros$Next_State=="0000")



states.list<-colnames(ast.twotab)

ast.twotab.m<-as.matrix(ast.twotab)
ast.MC<-new("markovchain",states=states.list,byrow=T,transitionMatrix=ast.twotab.m,name="Baseball States - Astros")




ast.exp.runs<-aggregate(Astros$Runs_Scored,by=list(Astros$State,Astros$Next_State),FUN=mean)
names(ast.exp.runs)<-c("State","Next_State","Exp_Runs")

tempmerge<-merge(ast.exp.runs,exp.runs,by=c("State","Next_State"),all=T)
for (m in 1:nrow(tempmerge)){
  if (is.na(tempmerge$Exp_Runs.x[m])){
    tempmerge$Exp_Runs.x[m]<-tempmerge$Exp_Runs.y[m]
  }
}
tempmerge<-tempmerge[,1:3]
ast.exp.runs<-tempmerge
names(ast.exp.runs)<-c("State","Next_State","Exp_Runs")

#Loop all possible states
all.states<-c("0000","0001","0010","0011","0100","0101","0110","0111","1000","1001","1010","1011","1100","1101","1110","1111","2000","2001","2010","2011","2100","2101","2110","2111","3000")
ast.exp.table<-NA
for (x in 1:25){
  i<-all.states[x]
  tot.sim.runs<-0
  for (c in 1:10000){
    i<-all.states[x]
    vec.states<-c(i)
    iter<-0
    current.exp.runs<-0
    while(i!="3000"){
      
      mci<-as.data.frame(ast.MC[i])
      names(mci)<-c("prob")
      mci["min"]<-0
      mci["max"]<-0
      mci["choice"]<-0
      mci$max[1]<-mci$prob[1]
      for (j in 2:nrow(mci)){
        k<-j-1
        mci$min[j]<-mci$max[k]
        mci$max[j]<-mci$min[j]+mci$prob[j]
      }
      mci$min<-mci$min*1000
      mci$max<-mci$max*1000
      #chance<-sample(1:1000,1)
      chance<-as.numeric(runif(1,min=0,max = 1000))
      for (a in 1:nrow(mci)){
        if ((chance>mci$min[a]) & (chance<=mci$max[a])){
          mci$choice[a]<-1
          n_state<-rownames(mci[a,])
        }
      }
      st.exp.runs<-subset(ast.exp.runs,ast.exp.runs$State==i & ast.exp.runs$Next_State==n_state)
      current.exp.runs<-current.exp.runs+st.exp.runs$Exp_Runs
      vec.states<-c(vec.states,n_state)
      i<-n_state
      iter<-iter+1
    }
    tot.sim.runs<-tot.sim.runs+current.exp.runs
  }
  avg.exp.runs<-tot.sim.runs/10000
  #print(all.states[x])
  #print(tot.sim.runs)
  #print(avg.exp.runs)
  row.add<-c(all.states[x],avg.exp.runs)
  ast.exp.table<-rbind(ast.exp.table,row.add)
}

ast.exp.table.clean<-ast.exp.table[-1,]
names(ast.exp.table.clean)<-c("State","Exp_Runs")
ast.exp.table.clean<-as.data.frame(ast.exp.table.clean)

##setwd("~/Northwestern University/Summer 2018 Quarter/MSDS 456/Assignments/Assignment 3/play_files")
#write.csv(exp.table.clean, "MLB_Expected_Runs.csv")
#write.csv(dod.exp.table.clean,"Dodgers_Expected_Runs.csv")
#write.csv(ast.exp.table.clean, "Astros_Expected-Runs.csv")
#write.csv(Astros, "Astros.csv")
#write.csv(Dodgers,"Dodgers.csv")
#write.csv(fullset,"MLB.csv")
#READ FILES IF STARTING HERE
#exp.table.clean<-read.csv("MLB_Expected_Runs.csv")
#dod.exp.table.clean<-read.csv("Dodgers_Expected_Runs.csv")
#ast.exp.table.clean<-read.csv("Astros_Expected-Runs.csv")
#Astros<-read.csv("Astros.csv")
#Dodgers<-read.csv("Dodgers.csv")
#fullset<-read.csv("MLB.csv")
#comp.table<-cbind(exp.table.clean[,1],exp.table.clean[,2],dod.exp.table.clean[,2],ast.exp.table.clean[,2])


##########################################
##########################################
##########################################
#THIS SECTION BUILDS TRANSITION MATRICES FOR EACH INDIVIDUAL PLAYER: THEY WILL BE CALLED "two,altuj001" , "two.brega001", etc.

dod.lineup.names<-c("Chris Taylor","Corey Seager","Justin Turner","Cody Bellinger","Yasiel Puig","Joc Pederson","Logan Forsythe","Austin Barnes","Pitcher")
ast.lineup.names<-c("Jose Altuve","Alex Bregman","Carlos Correa","Marwin Gonzalez","Yulleski Gurriel","Brian McCann","Josh Reddick","George Springer","Pitcher")
ast.lineup<-c("altuj001","brega001","corrc001","gonzm002","gurry001","mccab002","reddj001","sprig001")
dod.lineup<-c("taylc001","seagc001","turnj001","bellc002","puigy001","pedej001","forsl001","barna001")
lineups<-c("altuj001","brega001","corrc001","gonzm002","gurry001","mccab002","reddj001","sprig001","taylc001","seagc001","turnj001","bellc002","puigy001","pedej001","forsl001","barna001")


for (i in 1:8){
  pname<-dod.lineup[i]
  assign(pname,subset(Dodgers,Dodgers$Batter==as.character(pname)))
}

for (i in 1:8){
  pname<-ast.lineup[i]
  assign(pname,subset(Astros,Astros$Batter==as.character(pname)))
}
#2,5,6
player.list<-list(altuj001,brega001,corrc001,gonzm002,gurry001,mccab002,reddj001,sprig001,taylc001,seagc001,turnj001,bellc002,puigy001,pedej001,forsl001,barna001)
for (j in 1:16){
temp<-as.data.frame(player.list[j])
sub<-subset(temp,select = c("State","Next_State"))
twotab<-table(sub)
if (ncol(twotab)==25){  
  twentyfive<-as.data.frame(c(100,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0))
  p<-25
}
if (ncol(twotab)==24){  
  twentyfive<-as.data.frame(c(100,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0))
  p<-24
}
if (ncol(twotab)==23){  
  twentyfive<-as.data.frame(c(100,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0))
  p<-23
}
if (ncol(twotab)==22){  
  twentyfive<-as.data.frame(c(100,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0))
  p<-22
}
if (ncol(twotab)==21){  
  twentyfive<-as.data.frame(c(100,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0))
  p<-21
}
if (ncol(twotab)==20){  
  twentyfive<-as.data.frame(c(100,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0))
  p<-20
}
if (ncol(twotab)==19){  
  twentyfive<-as.data.frame(c(100,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0))
  p<-19
}

names(twentyfive)<-c("3000")
twentyfive<-t(twentyfive)

twotab<-rbind(twotab,twentyfive)
twotab<-as.data.frame(twotab)


twotab<-twotab[,1:p]
st.col<-25-ncol(twotab)
st.row<-25-nrow(twotab)

while (nrow(twotab)!=25) {
  k<-25-nrow(twotab)
  x<-26-ncol(twotab)
  if (x==1){
    radd<-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
  }
  if (x==2){
    radd<-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
  } 
  if (x==3){
    radd<-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
  } 
  if (x==4){
    radd<-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
  }
  for (i in 1:k){
      twotab<-rbind(twotab,radd)
    }
  }


  k<-25-ncol(twotab)
  #x<-26-nrow(twotab)
  #if (x==1){
    cadd<-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
  #}
  #if (x==2){
  #  cadd<-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
  #} 
  #if (x==3){
  #  cadd<-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
  #} 
  #if (x==4){
  #  cadd<-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
  #}
  if (k!=0){
    for (i in 1:k){
      twotab<-cbind(twotab,cadd)
    }
  }  

miss.row<-all.states[!all.states %in%rownames(twotab)]
miss.col<-all.states[!all.states %in%colnames(twotab)]
e<-length(miss.row)
f<-length(miss.col)

if (e!=0){
  for ( b in 1:e){
    row.names(twotab)[26-b]<-miss.row[b]
  }
}

if (f!=0){
  for ( b in 1:f){
    colnames(twotab)[26-b]<-miss.col[b]
  }
}

twotab<-twotab[c("0000","0001","0010","0011","0100","0101","0110","0111","1000","1001","1010","1011","1100","1101","1110","1111","2000","2001","2010","2011","2100","2101","2110","2111","3000")]
twotab<-twotab[c("0000","0001","0010","0011","0100","0101","0110","0111","1000","1001","1010","1011","1100","1101","1110","1111","2000","2001","2010","2011","2100","2101","2110","2111","3000"),]

twotab["Sum"]<-0
for (i in 1:25){
  twotab$Sum<-twotab$Sum+twotab[,i]
}

for (w in 1:25){
  if (twotab$Sum[w]<3){
    twotab[w,]<-mlb.twotab[w,]
  }
}

for (i in 1:25){
  twotab[,i]<-twotab[,i]/twotab$Sum
}
twotabx<-twotab[,1:25]
twotab<-as.matrix(twotabx)
twotab.temp<-twotab
twotab<-new("markovchain",states=states.list,byrow=T,transitionMatrix=twotab,name="Baseball States")

assign(paste("two.",lineups[j],sep = ""),twotab)

}


#MAKES A TRANSITION MATRIX FOR THE PITCHERS WHO WILL ALWAYS STRIKE OUT
two.pitcher<-twotab.temp
two.pitcher["0000",]<-c(0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
two.pitcher["0001",]<-c(0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
two.pitcher["0010",]<-c(0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
two.pitcher["0011",]<-c(0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0)
two.pitcher["0100",]<-c(0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0)
two.pitcher["0101",]<-c(0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0)
two.pitcher["0110",]<-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0)
two.pitcher["0111",]<-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0)
two.pitcher["1000",]<-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0)
two.pitcher["1001",]<-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0)
two.pitcher["1010",]<-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0)
two.pitcher["1011",]<-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0)
two.pitcher["1100",]<-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0)
two.pitcher["1101",]<-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0)
two.pitcher["1110",]<-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0)
two.pitcher["1111",]<-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0)
two.pitcher["2000",]<-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1)
two.pitcher["2001",]<-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1)
two.pitcher["2010",]<-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1)
two.pitcher["2011",]<-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1)
two.pitcher["2100",]<-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1)
two.pitcher["2101",]<-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1)
two.pitcher["2110",]<-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1)
two.pitcher["2111",]<-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1)
two.pitcher["3000",]<-c(1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)

two.pitcher<-new("markovchain",states=states.list,byrow=T,transitionMatrix=twotab.m,name="Baseball States")


#THIS IS WHERE THE GAME 7 SIMULATION BEGINS. I'VE SIMULATED TEN THOUSAND GAMES.
ast.lineup.twotab<-c(two.altuj001,two.brega001,two.corrc001,two.gonzm002,two.gurry001,two.mccab002,two.reddj001,two.sprig001,two.pitcher)
dod.lineup.twotab<-c(two.taylc001,two.seagc001,two.turnj001,two.bellc002,two.puigy001,two.pedej001,two.forsl001,two.barna001,two.pitcher)

tot.dod.score<-0
tot.ast.score<-0
ast.wins<-0
dod.wins<-0
for (n in 1:10000){
inning<-1
ast.score<-0
dod.score<-0
ast.spot<-1
dod.spot<-1
ast.current.exp.runs<-0
dod.current.exp.runs<-0
while (inning<10){
  #ASTROS AT-BAT
  cur.state<-"0000"
  while (cur.state!="3000"){
    batter<-ast.lineup.twotab[ast.spot][[1]]
    names(batter)<-c("0000","0001","0010","0011","0100","0101","0110","0111","1000","1001","1010","1011","1100","1101","1110","1111","2000","2001","2010","2011","2100","2101","2110","2111","3000")
    mci<-as.data.frame(batter[cur.state])
    names(mci)<-c("prob")
    mci["min"]<-0
    mci["max"]<-0
    mci["choice"]<-0
    mci$max[1]<-mci$prob[1]
    for (j in 2:nrow(mci)){
      k<-j-1
      mci$min[j]<-mci$max[k]
      mci$max[j]<-mci$min[j]+mci$prob[j]
    }
    mci$min<-mci$min*1000
    mci$max<-mci$max*1000
    #chance<-sample(1:1000,1)
    chance<-as.numeric(runif(1,min=0,max = 1000))
    for (a in 1:nrow(mci)){
      if ((chance>mci$min[a]) & (chance<=mci$max[a])){
        mci$choice[a]<-1
        n_state<-rownames(mci[a,])
      }
     }
 
  st.exp.runs<-subset(ast.exp.runs,ast.exp.runs$State==cur.state & ast.exp.runs$Next_State==n_state)
  ast.current.exp.runs<-ast.current.exp.runs+st.exp.runs$Exp_Runs
  cur.state<-n_state 
  #print(ast.spot)
  ast.spot<-ast.spot+1
  if (ast.spot==10){
    ast.spot<-1
  }
  }
  
###DODGERS AT-BAT
cur.state<-"0000"
while (cur.state!="3000"){
  batter<-dod.lineup.twotab[ast.spot][[1]]
  names(batter)<-c("0000","0001","0010","0011","0100","0101","0110","0111","1000","1001","1010","1011","1100","1101","1110","1111","2000","2001","2010","2011","2100","2101","2110","2111","3000")
  mci<-as.data.frame(batter[cur.state])
  names(mci)<-c("prob")
  mci["min"]<-0
  mci["max"]<-0
  mci["choice"]<-0
  mci$max[1]<-mci$prob[1]
  for (j in 2:nrow(mci)){
    k<-j-1
    mci$min[j]<-mci$max[k]
    mci$max[j]<-mci$min[j]+mci$prob[j]
  }
  mci$min<-mci$min*1000
  mci$max<-mci$max*1000
  #chance<-sample(1:1000,1)
  chance<-as.numeric(runif(1,min=0,max = 1000))
  for (a in 1:nrow(mci)){
    if ((chance>mci$min[a]) & (chance<=mci$max[a])){
      mci$choice[a]<-1
      n_state<-rownames(mci[a,])
    }
  }


st.exp.runs<-subset(dod.exp.runs,dod.exp.runs$State==cur.state & dod.exp.runs$Next_State==n_state)
dod.current.exp.runs<-dod.current.exp.runs+st.exp.runs$Exp_Runs
cur.state<-n_state 
dod.spot<-dod.spot+1
if (dod.spot==10){
  dod.spot<-1
}
}  
inning<-inning+1
#print(inning)
ast.score<-ast.current.exp.runs
dod.score<-dod.current.exp.runs
}


#print(ast.score)
#print(dod.score)
tot.ast.score<-tot.ast.score+ast.score
tot.dod.score<-tot.dod.score+dod.score

if (ast.score>dod.score){
  ast.wins<-ast.wins+1
  #print("ASTROS")
}

if (dod.score>ast.score){
  dod.wins<-dod.wins+1
  #print("DODGERS")
}
}


print(paste("Astros Average Runs Scored: ",tot.ast.score/10000))
print(paste("Dodgers Average Runs Scored: ",tot.dod.score/10000))
print(paste("Dodgers Total Wins: ",dod.wins))
print(paste("Astros Total Wins: ",ast.wins))


