classdef Logger < handle
    %LOGGER Stores log data! Simple as that
    %   Stores things in a huge matrix. So, must be storing columns. 
    
    properties
        names = {};
        storage = {};
        h = 0.1;
        time = 1;
        
    end
    
    methods
        
        function this = Logger()
           this = this@handle;
        end
        
        
    end
    
    methods
       
        function init(obj, time, h)
           obj.time = time;
           obj.h = h;
        end
        
        function add(obj, name, size)
           if nargin < 3
               error('Not enough!');
           end
           % Check if name is already present
           if isvalidname(obj, name)
               obj.names
               error(['Name ' name ' already present in logger!']);
           end
           obj.names{end+1} = name;
           obj.storage{end+1} = zeros(size, length(obj.time));
        end
        
        function set(obj, name, data)
           if ~isvalidname(obj, name)
                error('Invalid name');
           else

                obj.storage{obj.getId(name)} = data;
           end 
        end
        
        function store(obj, name, data, k)
           if ~isvalidname(obj, name)
                error('Invalid name');
           else
                obj.storage{obj.getId(name)}(:,k) = data;
           end 
        end
        
        % NB: Will copy the matrix
        function data = get(obj, name, index, time)
            
            if ~isvalidname(obj, name)
                error('Invalid name');
            else
                if nargin > 2
                    paren = @(x, varargin) x(varargin{:});
                    if nargin > 3
                        data = paren(obj.storage{getId(obj, name)}, index, time);
                    else
                        data = paren(obj.storage{getId(obj, name)}, index, ':');
                    end
                    
                        
                else                    
                   data = obj.storage{getId(obj, name)};
                end
            end
        end
        
        
        
        function ts = getts(obj, name)
            ts = timeseries(obj.get(name), obj.time, 'Name', name);
        end  
        
        function id = getId(obj, name)
            % Find ID from a name
            id = find(strcmp(obj.names, name));
        end
        
        function isvalid = isvalidname(obj,name)
            % checks if ID exists. 
            isvalid = ~isempty(getId(obj, name));
        end
        
        function stats = getInfo(obj)
           % Pulls all name and size info 
           stats = repmat(struct('name', [], 'size', 0), length(obj.names), 1);
           
           for i = 1:length(obj.names)
              stats(i).name =  obj.names{i};
              stats(i).size =  size(obj.storage{i},1); 
           end
           
        end
        
    end
    
end

