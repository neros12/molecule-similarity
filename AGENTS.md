# AGENTS.md

## Language Rules

- Use Korean for all user-facing natural-language responses.
- Preserve source code, identifiers, commands, file paths, names, quoted error messages, and logs in their original form.

## General Rules

- Prefer the smallest changes necessary to complete the task.
- Do not introduce abstractions for single-use code.
- Do not "improve" adjacent code, comments, formatting, or unrelated implementation.
- Do not search, inspect, or modify `.venv/`, `site-packages/`, or other installed dependency directories.

## Python Coding Rules

### General

- Use comprehensions only when the transformation for each item requires at most one simple operation or function call.
- Prefer human-readable formats such as `JSON` and `CSV` over binary formats when practical.
- Avoid `while` loops. Use bounded `for...else` loops, iterators, or explicit control-flow abstractions instead.
- For multi-line strings, always insert a newline immediately after the opening triple quotes (`'''` or `"""`).

## Blank Lines

### Before a Return After Processing

- Insert exactly one blank line before `return` when it follows calculations, assignments, or other processing.
- Apply this rule inside functions, conditional blocks, and loops, except for the immediate returns described below.

### Immediate Returns

- Do not insert a blank line between a function header and a `return` that is the first statement in its body.
- Do not insert a blank line between a control-flow block header and an immediate `return`.
- Calculations within the `return` expression do not count as preceding processing.

### Changes in Context

- Insert exactly one blank line when the purpose or processing stage changes.
- After an outermost `for` or `if` statement, insert exactly one blank line before the next statement at the same indentation level, except between consecutive `for` statements that serve the same purpose and processing stage.
- An outermost statement is not nested inside another `for` or `if` statement. Treat any attached `elif` or `else` clauses as part of the same statement, without inserting blank lines between those clauses.
- Keep consecutive calculations or settings that serve the same purpose together.
- A calculation and a simple check of its result may remain in the same group.
- Except for the outermost-statement rule above, do not add blank lines solely because a statement is an assignment, conditional, or loop.

## Type Safety

- Add type hints to function parameters.
- Do not explicitly annotate return types; let return types be inferred or propagated from the underlying function.
- Perform type conversions only when they are logically necessary and safe.
- Use `typing.cast()` only when the type is known from program logic but cannot be inferred by the type checker.
- Limit unavoidable `Any` and `object` types to untyped external-library or serialization boundaries.
- Do not introduce unresolved or `unknown` types.
- Before completing the task, run Pyright in `standard` mode and verify that no newly introduced type errors remain.

## Docstrings

- If writing docstring is requested:
  - Add docstrings to public classes, functions, and properties.
  - Write Python docstrings using the Google Python Style Guide format.
  - Begin each docstring with a concise summary written as an imperative sentence.
  - Include `Args`, `Returns`, `Raises`, and `Attributes` only when applicable.
  - Do not repeat information already clear from names or type annotations.
  - Keep section headings such as `Args`, `Returns`, `Raises`, and `Attributes` in English.

## Task Completion

Before completing a task:

- Verify that only necessary files and code were changed.
- Verify that unrelated code, comments, and formatting were left untouched.
- Run Pyright in `standard` mode and confirm that no newly introduced type errors remain for Python-related changes.
- After each Pyright check, delete any logs, reports, or temporary output files generated for that check.
- Do not run pytest or Ruff, including Ruff linting and formatting commands.

In the final response, include:

1. **Change description** — what was changed.
2. **Verification results** — tests, type checks, linting, or other validation performed.
3. **Follow-up notes** — remaining limitations, skipped verification, or relevant next steps.
