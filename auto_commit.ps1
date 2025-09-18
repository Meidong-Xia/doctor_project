# 进入仓库
Set-Location -Path $PSScriptRoot

# 检查是否有改动
$changes = git status --porcelain

if (-not [string]::IsNullOrWhiteSpace($changes)) {
    # 添加所有改动
    git add .

    # 提交（带时间戳）
    $time = Get-Date -Format "yyyy-MM-dd HH:mm:ss"
    git commit -m "Auto commit on $time"

    # 推送
    git push origin main   # 如果分支不是 main，改成你的分支名
} else {
    Write-Output "[$(Get-Date -Format "yyyy-MM-dd HH:mm:ss")] No changes to commit."
}
