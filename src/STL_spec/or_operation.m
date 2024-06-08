function  string  = or_operation(string1,string2)

%%% please consider the string '[i,j]' means the i-th predicate is imposed at time j an is valid. or in another word b_{phi_i}(s_j)>0. take a look at manual. 
%%%%%%%%% stacks the cell of char arrays and applies '|' between the two.

[size11,size21] = size(string1);
[size12,size22] = size(string2);
size1 = size11+size12;
size2 = max(size21,size22)+2;
str = cell(size1,size2);
str(:) = {' '};



if size21>=size22
    str(1,1) = {'('};
    str(1:size11, 2:size2-1) = string1;
    str(size11, size2-1 ) =  {')'};
    str(size11, size2   ) =  {'|'};
    shift =  0.5*(size2-size22);
    str(size11+1,shift) = {'('};
    str(size11+1:end, shift+1:size2-shift) = string2;
    str(size1, size2-shift ) =  {')'};
else
    str(1,1) = {'('};
    str(1:size12, 2:size2-1) = string2;
    str(size12, size2-1 ) =  {')'};
    str(size12, size2   ) =  {'|'};
    shift =  0.5*(size2-size21);
    str(size12+1,shift) = {'('};
    str(size12+1:end, shift+1:size2-shift) = string1;
    str(size1, size2-shift ) =  {')'};
end

string = str;

end