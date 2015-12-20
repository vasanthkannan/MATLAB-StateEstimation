% add_noise
% 
if numVmeas~=0
   for i=1:numVmeas
       if Vmeasstatus(i) == 1    
           if Vmeasbus(i) ~= refbus % do not assign noise to refbus
               Vmeasvaluenormrandom = Vmeassigma(i)*randn;
               Vmeasvalue(i) = Vmeasvalue(i)+Vmeasvaluenormrandom;
           end
       end
   end    
end

Slack = 1;
if numAmeas~=0
   for i=1:numAmeas
     if Ameasstatus(i) == 1
        if Ameasbus(i) ~= refbus  % do not assign noise to refbus
           Ameasvaluenormrandom = Ameassigma(i)*randn;
           Ameasvalue(i) = Ameasvalue(i)+Ameasvaluenormrandom;
        end
     end
   end    
end

if numImeas~=0
   for i=1:numImeas
       if Imeasstatus(i) == 1
           ImeasPvaluenormrandom = ImeasPsigma(i)*randn;
           ImeasPvalue(i) = ImeasPvalue(i)+ImeasPvaluenormrandom;
           ImeasQvaluenormrandom = ImeasQsigma(i)*randn;
           ImeasQvalue(i )= ImeasQvalue(i)+ImeasQvaluenormrandom;
       end
   end    
end
if numFmeas~=0
    for i=1:numFmeas
        if Fmeasstatus(i) == 1
            FmeasPvaluenormrandom=FmeasPsigma(i)*randn;
            FmeasPvalue(i)=FmeasPvalue(i)+FmeasPvaluenormrandom;
            FmeasQvaluenormrandom=FmeasQsigma(i)*randn;
            FmeasQvalue(i)=FmeasQvalue(i)+FmeasQvaluenormrandom;
        end
   end   
end
   
Noised=1;

fprintf('%s\n','noise added!');

Scadamenu
