% return back populations
function Moth_pos = Bound_Checking(Moth_pos,ub,lb)
for i=1:size(Moth_pos,1)
    Flag4ub=Moth_pos(i,:)>ub; 
    Flag4lb=Moth_pos(i,:)<lb;
    Moth_pos(i,:)=(Moth_pos(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
end 
        