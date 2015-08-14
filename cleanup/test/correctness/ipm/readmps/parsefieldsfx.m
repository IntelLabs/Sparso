%
% fields=parsefieldsfx(line)
%
% Returns a cell array of the fields in the input line, using MPS
% fixed format.  Assumes that the line coming in is a non-comment line.
%
function fields=parsefieldsfx(line)
%
% Setup a cell array of fields.
%
fields=cell(0);
%
% See if the first column is a space or not.  If not a space then
% this is a section header.  
%
if (line(1)==' ')
%
% It's a regular fixed format record, possibly with type in columns
% 2-3, followed by up to 3 fields, or without a type in columns
% 2-3, followed by up to 5 fields.
%
% First, figure out if there's a type in columns 2-3.
%
  if (strcmp(line(2:3),'  ')==1)
%
% Up to 5 fields in this line.
%
    if (length(line) >= 5)
      fields{1}=addspaces(line(5:min(12,length(line))),8);
    end
    if (length(line) >= 15)
      fields{2}=addspaces(line(15:min(22,length(line))),8);
    end
    if (length(line) >= 25)
      fields{3}=addspaces(line(25:min(36,length(line))),12);
    end
    if (length(line) >= 40)
      fields{4}=addspaces(line(40:min(47,length(line))),8);
    end
    if (length(line) >= 50)
      fields{5}=addspaces(line(50:min(61,length(line))),12);
    end
    if (length(line) >= 62)
      warning('Ignoring excess characters in line');
      line
    end
  else % not a type in columns 2-3.
%
% There is a type in columns 2-3.
%
    if (line(3)==' ')
      fields{1}=line(2:2);
    else
      if (line(2)==' ')
	fields{1}=line(3:3);
      else
	fields{1}=line(2:3);
      end
    end
    fields{1}=upper(fields{1});
%
% Now, get the remaining fields.
%
    if (length(line) >= 5)
      fields{2}=addspaces(line(5:min(12,length(line))),8);
    end
    if (length(line) >= 15)
      fields{3}=addspaces(line(15:min(22,length(line))),8);
    end
    if (length(line) >= 25)
      fields{4}=addspaces(line(25:min(36,length(line))),12);
    end
    if (length(line) >= 37)
      warning('Ignoring excess characters in line');
      line
    end

  end % line(2:3)='  '
else
%
% A section header.
%
  ptr=2;
  while ((line(ptr)~=' ') & (ptr < length(line)))
    ptr=ptr+1;
  end
  if (line(ptr)==' ')
    ptr=ptr-1;
  end
  fields{1}=line(1:ptr);

  if (length(line) >= 4)
    if ((strcmp(upper(line(1:4)),'NAME')==1) & (length(line)>=5))
      fields{2}=line(5:(min(12,length(line))));
    end
  end
end
%
% Ignore any trailing fields of all blanks.
%
for i=length(fields):-1:1
  if (strcmp(fields{i},'        ')==1)
    fields{i}=[];
    continue
  end
  if (strcmp(fields{i},'            ')==1)
    fields{i}=[];
    continue
  end
  break
end
%
% Get rid of any empty fields.
%
for i=length(fields):-1:1
  if (~isempty(fields{i}))
    l=i;
    break
  end
end
fields=fields(1:l);

