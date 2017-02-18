function Z=mpoweroff(X,Y)
%POWER Summary of this function goes here
%   Detailed explanation goes here    
    Z=realpow(X,Y);
%     msgID='mpower:complexInput';
%     msgText='(%.4g)^(%.4g) is complex.\n>\tIn %s (line ';
%     try
%         if ~isreal(builtin('mpower',X,Y))
%             ME=MException(msgID,msgText);
%             throw(ME);
%         else
%             Z=builtin('mpower',X,Y);
%         end
%     catch ME
%         %warning('off','backtrace')
%         fileroot=ME.stack(end).name;
%         file=ME.stack(2).name;
%         filepath=[fileroot '.m'];
%         line=ME.stack(2).line;
%         %warning(msgText,X,Y,file,line,...
%            % fprintf('<a href="matlab:opentoline(%s,%d,0)">GoTo Line %d</a>',file,line,line));
%         str='(%.4g)^(%.4g) is complex.\n>\tIn %s (line <a href="matlab:opentoline(%s,%d,0)">%d</a>)\n';
%         warning(str,X,Y,file,filepath,line,line);
%         Z=NaN;
%     end

end

