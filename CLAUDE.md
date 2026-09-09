# biophys_plotting_agent

This repo IS the source of the `biophys-plotting` Claude Code plugin (marketplace `biophys-tools`,
installed from `github.com/tru489/biophys_plotting_agent`). Slash commands and skills invoked as
`biophys-plotting:*` run from a separate cached copy in `~/.claude/plugins/cache/`, not from this
working directory — editing files here has no effect on running sessions until that cache is
refreshed.

## After completing any change to this plugin

Whenever a change here is finished (code, skill, or doc), update the local install so future
slash/skill invocations use the current version:

1. Bump `version` in `.claude-plugin/plugin.json` (semver — patch for fixes, minor for features).
2. Commit, ending the commit message with `version: X.Y.Z` (matches this repo's existing commit
   history convention — see `git log`).
3. Push to `origin main` (the marketplace source is this GitHub repo).
4. Run `claude plugin marketplace update biophys-tools` then
   `claude plugin update biophys-plotting@biophys-tools` to refresh the local cache.
5. Tell the user a restart of Claude Code is needed for the update to take effect.

Do this automatically at the end of a change, without waiting to be asked each time.
