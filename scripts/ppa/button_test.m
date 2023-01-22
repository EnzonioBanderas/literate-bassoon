function [outputArg1,outputArg2] = untitled3(inputArg1,inputArg2)
ax_screen = subplot(5,2,1:4)

button_prev = uicontrol("Style","pushbutton","String","Previous","Position",[0 0 60 20])
button_prev.Callback = @prevButtonPush
button_continue = uicontrol("Style","pushbutton","String","Next","Position",[60 0 60 20])
button_stop = uicontrol("Style","pushbutton","String","Stop","Position",[120 0 60 20])

x = 2;
function prevButtonPush(~,~)
    x = x-1
end

end