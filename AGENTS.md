<!-- Carter writes this file once and does not overwrite it. -->

# dealii-rust

## Workspace Configuration

- **Mode**: umbrella

## Agent System

This workspace uses the Carter multi-agent system. Requests are routed to specialist agents based on domain.

### Available Agents

- `carter-maestro-local`
- `carter-architect-local`
- `carter-engineer-local`
- `carter-devops-local`
- `carter-pm-local`
- `carter-learner-local`
- `carter-narrator-local`
- `carter-security-local`

### Routing

Carter Maestro (`carter-maestro-local`) is the orchestrator. Direct specialist requests to the appropriate agent. When uncertain, Maestro routes.

## Conventions

- **Knowledge**: two sets, read before acting — workspace first, then global.
  - `.carter/knowledge/` — this workspace's files (gotchas, patterns, decisions).
  - `~/.carter/knowledge/` — global seeds installed by `carter install` and grown by `carter-distill`; shared across workspaces.
- **State**: `.carter/state/` contains plans and session checkpoints.
- **Plans**: Active plans in `.carter/state/plans/active/`. Check at session start.
- **Current date**: Obtain the current date from a tool or the environment (a
  native OpenCode date variable where available, otherwise run `date`)
  **before writing any dated frontmatter** (`created:`/`decided:`/`updated:`). A
  training-cutoff guess is wrong; never bake or guess a date.
- **TDD**: Spec-handoff convention enabled. Document contracts before complex work.
- **Quality Gate**: Ask the maestro in chat (e.g. "review", "check before push"). Max 2 iterations.
- **Redaction**: External-facing content goes through carter-narrator-local for redaction.
- **Language**: Repository and outward content (commits, branches, code comments, docs, PRs, issues, Jira, Slack) is always English. Carter internal docs (ADRs, RFCs, `.carter/`, memories) use `en` (`docs_language`). Talk to the user in their language.
- **Commits**: The user is the author and committer; never pass `--author`, `--no-verify`, `-s` or `-S`. Credit comes only from the `coauthor_trailer` that `carter ops config --emit --field coauthor_trailer` prints (use the `carter-ops` skill); if `commits` is `refused`, stop and report.
- **Pushing**: Pushing is human-only: tell the user: Run `carter push` yourself, in your terminal or as `!carter push` in OpenCode, not through an agent. Only carter-devops-local may ask with `carter_push`; the user approves each use. After the push, `gh pr create` is fine.

## Security

- Always consider OWASP Top 10 in code and architecture reviews.
- Never hardcode secrets. Flag insecure patterns immediately.
- Validate inputs at system boundaries.

## Implementation Discipline

- After every change (feature, fix, refactor), remove what is no longer used — dead code, unused variables, orphaned helpers, stale imports.
- Reuse existing abstractions (classes, variables, functions, methods, interfaces) rather than reimplementing or hardcoding values already defined elsewhere.
- Do this cleanup as part of the same change, not as a follow-up.
- Before finishing, grep for each removed/replaced reference to confirm nothing else depends on it.

## Workflow

- Propose tests before implementation unless explicitly told otherwise.
- Think before acting — read the existing codebase and documentation first.
- Accuracy is more important than brevity, but avoid unnecessary verbosity.

## Technology Defaults

- Languages: Python, Go, bash
- IaC: Terraform
- CI/CD: GitHub Actions
- Python tooling: uv
- Version control: Git, GitHub
