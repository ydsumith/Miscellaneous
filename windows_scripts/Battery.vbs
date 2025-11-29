Option Explicit

'========================================
' Declarations & constants
'========================================
Dim objWMI, colProcesses, objProcess
Dim thisScript, procCount, answer, maxProcId

Dim oLocator, oServices, oResults, oResult
Dim iFull, iRemaining, bCharging, iPercent
Dim alreadyFullNotified
Dim isNearlyFull
Dim statusText, actionText, powerText

Const FULL_PERCENT_THRESHOLD = 98      ' % at which we consider "full"
Const FULL_RATIO_THRESHOLD   = 0.98    ' 98% of full capacity
Const LOW_BATTERY_THRESHOLD  = 20      ' low battery threshold (%)
Const POLL_INTERVAL_MS       = 60000   ' 60 seconds

alreadyFullNotified = 0

'========================================
' Single-instance protection
' - Keep newest instance, optionally kill older ones
'========================================
thisScript = LCase(WScript.ScriptFullName)
procCount = 0
maxProcId = 0

Set objWMI = GetObject("winmgmts:\\.\root\cimv2")
Set colProcesses = objWMI.ExecQuery( _
    "SELECT * FROM Win32_Process WHERE Name='wscript.exe' OR Name='cscript.exe'")

For Each objProcess In colProcesses
    On Error Resume Next
    Dim cmd
    cmd = LCase(objProcess.CommandLine)
    On Error GoTo 0

    If cmd <> "" Then
        If InStr(1, cmd, thisScript, vbTextCompare) > 0 Then
            procCount = procCount + 1
            If objProcess.ProcessId > maxProcId Then
                maxProcId = objProcess.ProcessId ' assume newest has highest PID
            End If
        End If
    End If
Next

If procCount > 1 Then
    answer = MsgBox( _
        "Another instance of this Battery Monitor script is already running." & vbCrLf & vbCrLf & _
        "Do you want to CLOSE the existing instance(s) and continue with this one?" & vbCrLf & _
        "(Yes = close old and keep this one, No = exit this one)", _
        vbQuestion + vbYesNo + vbDefaultButton2, "Battery Monitor")

    If answer = vbYes Then
        ' Kill all other instances of this script, keep the newest (this one)
        For Each objProcess In colProcesses
            On Error Resume Next
            cmd = LCase(objProcess.CommandLine)
            On Error GoTo 0

            If cmd <> "" Then
                If InStr(1, cmd, thisScript, vbTextCompare) > 0 Then
                    ' Only terminate if it's NOT the newest instance
                    If objProcess.ProcessId <> maxProcId Then
                        On Error Resume Next
                        objProcess.Terminate
                        On Error GoTo 0
                    End If
                End If
            End If
        Next
    Else
        WScript.Quit
    End If
End If

'========================================
' Battery monitor setup
'========================================
Set oLocator = CreateObject("WbemScripting.SWbemLocator")
Set oServices = oLocator.ConnectServer(".", "root\wmi")

'-------------------------------
' Get full charged capacity
'-------------------------------
iFull = 0
Set oResults = oServices.ExecQuery("SELECT * FROM BatteryFullChargedCapacity")

For Each oResult In oResults
    If Not IsNull(oResult.FullChargedCapacity) Then
        iFull = oResult.FullChargedCapacity
        Exit For
    End If
Next

If iFull <= 0 Then
    MsgBox "Unable to read battery full-charged capacity. Exiting.", _
           vbCritical, "Battery Monitor"
    WScript.Quit
End If

'-------------------------------
' Initial battery status
'-------------------------------
Set oResults = oServices.ExecQuery("SELECT * FROM BatteryStatus")
For Each oResult In oResults
    iRemaining = oResult.RemainingCapacity
    bCharging  = oResult.Charging  ' True = charger connected
Next

If iFull > 0 Then
    iPercent = Int((iRemaining / iFull) * 100)
Else
    iPercent = 0
End If

If iPercent < 0 Then iPercent = 0
If iPercent > 100 Then iPercent = 100

' Determine if battery is nearly full at startup
isNearlyFull = False
If iPercent >= FULL_PERCENT_THRESHOLD Then
    isNearlyFull = True
ElseIf iRemaining >= (FULL_RATIO_THRESHOLD * iFull) Then
    isNearlyFull = True
End If

'========================================
' First-run status message for the user
'========================================
If bCharging Then
    powerText = "Power: Charger is connected (charging)."
Else
    powerText = "Power: Running on battery (not charging)."
End If

If bCharging Then
    If isNearlyFull Then
        actionText = "Action: Battery is almost full. You can remove the charger if you wish."
    Else
        actionText = "Action: No immediate action needed. Charging in progress."
    End If
Else
    If iPercent < LOW_BATTERY_THRESHOLD Then
        actionText = "Action: Battery is low. Please connect the charger."
    Else
        actionText = "Action: No immediate action needed. You are on battery power."
    End If
End If

statusText =  "Battery Monitor started." & vbCrLf & vbCrLf & _
              "Battery: " & iPercent & "%" & vbCrLf & _
              powerText & vbCrLf & _
              actionText & vbCrLf & vbCrLf & _
              "Plugin status: Monitoring battery every " & (POLL_INTERVAL_MS \ 1000) & " seconds."

MsgBox statusText, vbInformation, "Battery Monitor - Initial Status"

'========================================
' Main monitoring loop
'========================================
Do While True
    ' Refresh battery status
    Set oResults = oServices.ExecQuery("SELECT * FROM BatteryStatus")
    For Each oResult In oResults
        iRemaining = oResult.RemainingCapacity
        bCharging  = oResult.Charging
    Next

    If iFull > 0 Then
        iPercent = Int((iRemaining / iFull) * 100)
    Else
        iPercent = 0
    End If

    If iPercent < 0 Then iPercent = 0
    If iPercent > 100 Then iPercent = 100

    ' Nearly full check
    isNearlyFull = False
    If iPercent >= FULL_PERCENT_THRESHOLD Then
        isNearlyFull = True
    ElseIf iRemaining >= (FULL_RATIO_THRESHOLD * iFull) Then
        isNearlyFull = True
    End If

    '------------------------------------
    ' Full charge notification (only if charging)
    '------------------------------------
    If bCharging And isNearlyFull And alreadyFullNotified = 0 Then
        MsgBox "Battery is at " & iPercent & "%." & vbCrLf & _
               "You can safely remove the charger.", _
               vbInformation, "Battery Monitor"
        alreadyFullNotified = 1
    End If

    ' Reset flag if we drop a bit below "full"
    If iPercent < (FULL_PERCENT_THRESHOLD - 5) Then
        alreadyFullNotified = 0
    End If

    '------------------------------------
    ' Low battery warning (only if NOT charging)
    '------------------------------------
    If (Not bCharging) And iPercent < LOW_BATTERY_THRESHOLD Then
        MsgBox "Battery is low (" & iPercent & "%)." & vbCrLf & _
               "Please connect the charger!", _
               vbExclamation, "Battery Monitor"
    End If

    ' Wait before checking again
    WScript.Sleep POLL_INTERVAL_MS
Loop
