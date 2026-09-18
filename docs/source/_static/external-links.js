/* Open off-site links in a new tab.
 *
 * Sphinx tags them class="reference external", but that class is also given to
 * absolute links back into this same site, so the hostname is what is tested
 * here.  rel="noopener" is set with it: without it the opened page gets a
 * window.opener handle back to this one.
 */
document.addEventListener('DOMContentLoaded', function () {
    var here = window.location.hostname;
    var links = document.querySelectorAll('a[href^="http://"], a[href^="https://"]');
    Array.prototype.forEach.call(links, function (a) {
        if (a.hostname && a.hostname !== here) {
            a.setAttribute('target', '_blank');
            a.setAttribute('rel', 'noopener noreferrer');
        }
    });
});
