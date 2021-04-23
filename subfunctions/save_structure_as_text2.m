function save_structure_as_text2(s,filename)

%save a structure s into the file filename as a text file
%supports substructures (but not subsubstructures)

fid1 = fopen(filename,'w');
names = fieldnames(s);
vals = struct2cell(s);

for i=1:size(names,1)
    
    if ~isstruct(vals{i,1})
        if size(vals{i,1},1)>1
            if ~iscell(vals{i,1})
                str = [names{i,1},' (truncated): ',num2str(vals{i,1}(1,:))];
                fprintf(fid1,'%s \r\n',str);
            else
                for j=1:size(vals{i,1},1)
                    if j ==1
                        strv = [names{i,1},': ',vals{i,1}{j,1}];
                    else
                        strv = [strv; [names{i,1},': ',vals{i,1}{j,1}]];
                    end
                end
                strv = cellstr(strv);
                fprintf(fid1,'%s\r\n',strv{:});
                
            end
        else
            if isa(vals{i,1},'numeric')
                str = [names{i,1},': ',num2str(vals{i,1})];
            elseif isa(vals{i,1},'char')
                str = [names{i,1},': ',vals{i,1}];
            else
                str = [names{i,1},': could not read value'];
            end
            fprintf(fid1,'%s \r\n',str);
        end
        
         
    else
        s2 = vals{i,1};
        names2 = fieldnames(s2);
        vals2 = struct2cell(s2);
        for j = 1:size(names2,1)
            if size(vals2{j,1},1) >1
                if ~iscell(vals2{j,1})
                    str = [names{i,1},'.',names2{j,1},' truncated: ',num2str(vals2{j,1}(1,:))];
                    fprintf(fid1,'%s \r\n',str);
                else
                    fprintf(fid1,'%s\r\n%s\r\n',[names{i,1},'.',names2{j,1}],vals{i,1}{:});
                end
            else
                str = [names{i,1},'.',names2{j,1},': ',num2str(vals2{j,1})];
            end
            
            fprintf(fid1,'%s \r\n',str);
        end
    end
end
fclose(fid1); 
end
  


