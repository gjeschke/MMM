function webcall( url, option )
%function webcall( url, option )
%   frame for call of Matlab's web function to select between system and
%   Matlab browser

global browser

if nargin<2,
    option='-none';
end;

switch lower(browser)
    case 'system'
        web(url,'-browser');
    case 'help'
        web(url,'-helpbrowser');
    case 'matlab'
        web(url);
    case 'coded'
        if nargin>1,
            web(url,option);
        else
            web(url,'-notoolbar');
        end;
    case 'notoolbar'
        web( url);
    case 'mixed'
        if strcmp(option,'-helpbrowser'),
            web(url,'-helpbrowser');
        else
            web(url,'-browser');
        end;
    otherwise
        if nargin>1,
            web(url,option);
        else
            web(url);
        end;
end

