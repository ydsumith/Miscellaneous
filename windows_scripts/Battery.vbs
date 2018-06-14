set oLocator = CreateObject("WbemScripting.SWbemLocator")
set oServices = oLocator.ConnectServer(".","root\wmi")
set oResults = oServices.ExecQuery("select * from batteryfullchargedcapacity")
for each oResult in oResults
iFull = oResult.FullChargedCapacity
next

Dim already
already = 0

set oResults = oServices.ExecQuery("select * from batterystatus")
for each oResult in oResults
iRemaining = oResult.RemainingCapacity
bCharging = oResult.Charging
next
iPercent = ((iRemaining / iFull) * 100) mod 100

msgbox "Currently Battery is at " & iPercent & "%",vbInformation, "Battery monitor"

while already = 0
set oResults = oServices.ExecQuery("select * from batterystatus")
for each oResult in oResults
iRemaining = oResult.RemainingCapacity
bCharging = oResult.Charging
next
iPercent = ((iRemaining / iFull) * 100) mod 100
If bCharging and (iPercent >= 98) Then 
	If already = 0 Then
		msgbox "Battery is at " & iPercent & "%, Remove the Charger!",vbInformation, "Battery monitor"
		already = 1
	End If
End If
wscript.sleep 60000 ' 1 minutes = 1 x 60s x 1000ms
wend
