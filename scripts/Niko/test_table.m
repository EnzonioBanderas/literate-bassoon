T = table(rand(2,1),rand(2,1),'VariableNames',{'Age','Weight'},'RowNames',{'P1','P2'}) ;
for i = 3:10
    Tnew = table(rand,rand,'VariableNames',{'Age','Weight'},'RowNames',{['P',num2str(i)]}) ;
    T = [T ; Tnew] 
end