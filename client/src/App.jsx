import React from 'react';
import {
  BrowserRouter as Router,
  Link,
  Route,
  Routes,
  useLocation,
} from 'react-router-dom';
import Home from './Home';
import Allosteric from './Allosteric';
import PDBCompare from './PDBCompare';
import NewAllosteric from './NewAllosteric';
import CholNet from './CholNet';
import './App.css';

function App() {
  return (
    <Router>
      <AppContent />
    </Router>
  );
}

function AppContent() {
  const location = useLocation();
  const isHomePage = location.pathname === '/';

  return (
    <>
      <Routes>
        <Route path="/" element={<Home />} />
        <Route path="/allosteric" element={<Allosteric />} />
        <Route path="/pdbpaircompare" element={<PDBCompare />} />
        <Route path="/newallosteric" element={<NewAllosteric />} />
        <Route path="/cholbindnet" element={<CholNet />} />
      </Routes>

      {!isHomePage && <NavBar />}
    </>
  );
}

function NavBar() {
  const location = useLocation();

  const links = [
    { to: '/', label: 'Home' },
    { to: '/newallosteric', label: 'Allosteric' },
    { to: '/pdbpaircompare', label: 'PDB Compare' },
    { to: '/cholbindnet', label: 'CholBindNet' },
  ];

  return (
    <nav className="app-nav" aria-label="Tool navigation">
      {links
        .filter((link) => link.to !== location.pathname)
        .map((link) => (
          <Link key={link.to} to={link.to}>
            {link.label}
          </Link>
        ))}
    </nav>
  );
}

export default App;
          