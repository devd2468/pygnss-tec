# MATLAB Cleanup Script
# Run this with Administrator privileges
Write-Host "MATLAB Cleanup Script" -ForegroundColor Cyan
Write-Host "======================" -ForegroundColor Cyan

# Stop any running MATLAB processes
Write-Host "`nStopping MATLAB processes..." -ForegroundColor Yellow
Get-Process matlab -ErrorAction SilentlyContinue | Stop-Process -Force

# Uninstall MATLAB via Windows Installer
Write-Host "`nAttempting to uninstall MATLAB..." -ForegroundColor Yellow
$matlabPath = "C:\Program Files\MATLAB"

# Remove folders
Write-Host "Removing MATLAB directories..." -ForegroundColor Yellow
$folders = @(
    "$matlabPath",
    "$env:APPDATA\MATLAB",
    "$env:LOCALAPPDATA\MATLAB",
    "$env:USERPROFILE\Documents\MATLAB",
    "$env:APPDATA\MathWorks"
)

foreach ($folder in $folders) {
    if (Test-Path $folder) {
        try {
            Remove-Item -Path $folder -Force -Recurse -ErrorAction Stop
            Write-Host "Removed: $folder" -ForegroundColor Green
        } catch {
            Write-Host "Failed to remove: $folder" -ForegroundColor Red
        }
    } else {
        Write-Host "Not found: $folder" -ForegroundColor Gray
    }
}

# Clean registry entries
Write-Host "`nCleaning registry..." -ForegroundColor Yellow
$regPaths = @(
    "HKLM:\SOFTWARE\MATLAB",
    "HKLM:\SOFTWARE\WOW6432Node\MATLAB",
    "HKCU:\SOFTWARE\MATLAB",
    "HKCU:\SOFTWARE\MathWorks"
)

foreach ($regPath in $regPaths) {
    if (Test-Path $regPath) {
        try {
            Remove-Item -Path $regPath -Force -Recurse -ErrorAction Stop
            Write-Host "Removed: $regPath" -ForegroundColor Green
        } catch {
            Write-Host "Failed to remove: $regPath" -ForegroundColor Red
        }
    } else {
        Write-Host "Not found: $regPath" -ForegroundColor Gray
    }
}

# Clean PATH environment variables
Write-Host "`nCleaning PATH..." -ForegroundColor Yellow
$currentPath = [Environment]::GetEnvironmentVariable("PATH", "Machine")
if ($currentPath -match "MATLAB") {
    $newPath = $currentPath -replace ";*C:\\Program Files\\MATLAB[^;]*", ""
    [Environment]::SetEnvironmentVariable("PATH", $newPath, "Machine")
    Write-Host "Cleaned PATH" -ForegroundColor Green
}

Write-Host "`nCleanup complete!" -ForegroundColor Cyan
Write-Host "You can now install MATLAB R2025b" -ForegroundColor Green
