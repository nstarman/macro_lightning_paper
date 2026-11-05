# The published figures, copied as they are. showyourwork's own copy of src/static breaks on a
# path with a space; a rule here takes precedence over it.
rule static_figures:
    input:
        "src/static/macro_schematic.pdf",
        "src/static/lightningconstraints.pdf",
    output:
        "src/tex/figures/macro_schematic.pdf",
        "src/tex/figures/lightningconstraints.pdf",
    shell:
        "cp {input[0]:q} {output[0]:q} && cp {input[1]:q} {output[1]:q}"

# The paper's docs site (syw-site).
try:
    import syw_site
except ImportError:
    syw_site = None
    print("syw-site is not installed in showyourwork's environment: the docs rules are off.")
if syw_site is not None:
    include: syw_site.rules_path()
