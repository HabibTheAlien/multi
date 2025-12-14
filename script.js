document.addEventListener('DOMContentLoaded', () => {
    const downloadRLink = document.getElementById('downloadR');
    const downloadPDFLink = document.getElementById('downloadPDF');
    const messageElement = document.getElementById('message');

    // Function to handle the download click
    const handleDownload = (fileName) => {
        // Display a confirmation message
        messageElement.textContent = `Download started for ${fileName}! Thank you.`;
        
        // Optionally, clear the message after a few seconds
        setTimeout(() => {
            messageElement.textContent = '';
        }, 5000); 
    };

    // Add event listeners to both links
    downloadRLink.addEventListener('click', () => {
        handleDownload('stata.txt');
    });

    downloadPDFLink.addEventListener('click', () => {
        handleDownload('text.pdf');
    });
});