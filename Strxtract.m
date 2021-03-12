function [ out ] = Strxtract( instr,indx )


inic=find(instr=='[',indx);
endc=find(instr==']',indx);
out=instr(inic+1:endc-1);

end
