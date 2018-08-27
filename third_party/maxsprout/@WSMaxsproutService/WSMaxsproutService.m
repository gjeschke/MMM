function obj = WSMaxsproutService

obj.endpoint = 'http://www.ebi.ac.uk/Tools/es/ws-servers/WSMaxsprout';
obj.wsdl = 'http://www.ebi.ac.uk/Tools/webservices/wsdl/WSMaxsprout.wsdl';

obj = class(obj,'WSMaxsproutService');

