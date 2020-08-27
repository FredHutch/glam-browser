function onClick() {
    dataLayer.push({ 'event': 'interaction' });
}

window.addEventListener("click", function () { onClick() });
