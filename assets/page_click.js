function onClick() {
    dataLayer.push({ 'event': 'interaction' });
}

window.dataLayer = window.dataLayer || [];
window.addEventListener("click", function () { onClick() });
