/**
 * Boolean parameters that survive the Nextflow command line.
 *
 * WHY THIS CLASS EXISTS
 * Nextflow 26.04 passes `--flag false` through as the STRING "false", and every
 * non-empty String is truthy in Groovy. So `--skip_metabuli false` leaves the skip
 * ON while reading as though it turned it off - the flag can only ever be set TRUE
 * from the command line. A -params-file entry parses as a real boolean and behaves
 * correctly, so the two routes silently disagree.
 *
 * That is not hypothetical here. scripts/run_bakeoff.sh carried
 * `plus_krakenuniq|--skip_krakenuniq false`, written specifically to repair an
 * earlier matrix where the variant was a byte-for-byte duplicate of the baseline.
 * The replacement spelling did not work either, so decisions.md D2 - which names
 * that exact variant as the test that would reverse the Kraken2-primary decision -
 * stayed unfalsifiable by its own stated test, twice over.
 *
 * params is read-only at runtime (assigning params[k] is silently ignored), so the
 * coercion cannot be done in place. Every read of a boolean parameter goes through
 * here instead.
 */
class BtfBool {

    private static final List<String> TRUE_WORDS  = ['true',  'yes', 'y', '1', 'on']
    private static final List<String> FALSE_WORDS = ['false', 'no',  'n', '0', 'off', '']

    /** Coerce a parameter value to a boolean, or throw if it is not boolean-ish. */
    static boolean of(Object v) {
        if (v == null)             return false
        if (v instanceof Boolean)  return (Boolean) v
        if (v instanceof Number)   return ((Number) v).intValue() != 0
        if (v instanceof CharSequence) {
            String t = v.toString().trim().toLowerCase()
            if (t in TRUE_WORDS)  return true
            if (t in FALSE_WORDS) return false
            throw new IllegalArgumentException(
                "not a boolean: ${v.inspect()} - use true/false, or omit the flag for the default")
        }
        throw new IllegalArgumentException(
            "not a boolean: ${v.inspect()} (${v.getClass().simpleName})")
    }

    /** Same, but names the parameter in the error so the message is actionable. */
    static boolean of(String name, Object v) {
        try { return of(v) }
        catch (IllegalArgumentException e) { throw new IllegalArgumentException("--${name}: ${e.message}") }
    }
}
