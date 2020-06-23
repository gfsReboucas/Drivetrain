classdef SimpackCOM
    %SIMPACKCOM wrapper
    properties(Access = private)
        COM;
    end
    
    properties
        mode; % GUI, SoLVer, Post
    end
    
    methods
        function obj = SimpackCOM(varargin)
            default = {'mode', 'SLV', ...
                       'version', ''};
            [flag, msg] = SimpackCOM.is_installed();
            default = process_varargin(default, varargin);
            
            if(flag)
                obj.COM = actxserver(sprintf('Simpack.%s%s', default.mode, ...
                                                      default.version));
            else
                error(msg.identifier, "%s", msg.message);
            end
            
            obj.mode = default.mode;
            
        end
    end
    
    methods(Static)
        function [val, msg] = is_installed()
            val = true;
            msg = [];
            try
                actxserver('Simpack.Post');
            catch msg
                val = false;
            end
        end
        
        function view_COM_methods()
            [flag, msg] = KISSsoftCOM.is_installed();
            if(flag)
                methodsview(actxserver('Simpack.Gui'));
                methodsview(actxserver('Simpack.Slv'));
                methodsview(actxserver('Simpack.Post'));
            else
                warning(msg.identifier, "%s", msg.message);
            end
        end
        
    end
    
end
